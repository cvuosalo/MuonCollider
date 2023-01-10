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

  fxHistName = GetString("xHistName", "x");
  fPositionHistName = GetString("PositionHistName", "z_r");
  fPhiHistName = GetString("PhiHistName", "phi");
  fThetaHistName = GetString("ThetaHistName", "theta");
  fPdgEnergyHistName = GetString("PdgEnergyHistName", "pdgid_energy");
  fMomentumHistName = GetString("MomentumHistName", "px_py_pz");
  
  fileName = GetString("FileName", "file");

  fInputArray = ImportArray(GetString("InputArray", "Delphes/allParticles"));
  fStableInputArray = ImportArray(GetString("StableInputArray", "Delphes/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();
  fItStableInputArray = fStableInputArray->MakeIterator();

  fOutputArray = ExportArray(GetString("OutputArray", "allParticles"));
  fStableOutputArray = ExportArray(GetString("StableOutputArray", "stableParticles"));

  file = TFile::Open(fileName);
  
  fxHist = (TH1F*) file->Get(fxHistName);
  fPosistionHist = (TH2F*) file->Get(fPositionHistName);
  fPhiHist = (TH1F*) file->Get(fPhiHistName);
  fThetaHist = (TH1F*) file->Get(fThetaHistName);
  fPdgEnergyHist = (TH2F*) file->Get(fPdgEnergyHistName);
  fMomentumHist = (TH3F*) file->Get(fMomentumHistName);
}

void BIBModule::Finish()
{
  if(fItInputArray) delete fItInputArray; 
  if(fItStableInputArray) delete fItStableInputArray; 
}

void Process()
{
  GenParticle *original, *bib;
  Float_t px, py, pz, energy, theta, eta, phi, m, charge, pdgid, x, y, z, r;

  if (!fxHist) return;
  if (!fPositionHist) return;
  if (!fPhiHist) return;
  if (!fThetaHist) return;
  if (!fMomentumHist) return;
  if (!fPdgEnergyHist) return;

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
    fPhiHist -> GetRandom(phi);
    fThetaHist -> GetRandom(theta);
    fxHist -> GetRandom(x);
    fPositionHist -> GetRandom2(z, r);
    fMomentumHist -> GetRandom3(px, py, pz);
    fPdgEnergyHist -> GetRandom2(pdgid, energy);
    
    if(energy <= 0.0) continue;
    
    bib.X = x;
    y = TMath::Sqrt(r * r - x * x);
    bib.Y = y;
    bib.Z = z;
    eta = - TMath::Log(TMath::Tan(theta/2));
    bib.Eta = eta;
    bib.Phi = phi;
    bib.E = energy;
    bib.Px = px;
    bib.Py = py;
    bib.Pz = pz;

    bib.PID = pdgid; 

    //assign mass and charge by pdgID here

    if (TMath::Abs(pdgid) == 11) {
        //electron
	if (pdgid > 0) {
	    charge = -1;
	} else {
	    charge = +1;
	}
	bib.Charge = charge;
	mass = 0.00051099895;
	bib.Mass = mass;
    } else if (TMath::Abs(pdgid) == 13) {
        //muon
	if (pdgid > 0) {
	    charge = -1;
	} else {
	    charge = +1;
	}
	charge = -1;
	bib.Charge = charge;
	mass = 0.1056583755;
	bib.Mass = mass;
    } else if (TMath::Abs(pdgid) == 22) {
	//photon
	charge = 0;
	bib.Charge = charge;
	mass = 0;
	bib.Mass = mass;
    } else if (TMath::Abs(pdgid) == 111) {
        //neutral pion
	charge = 0;
	bib.Charge = charge;
	mass = 0.1349768;
	bib.Mass = mass;
    } else if (TMath::Abs(pdgid) == 211) {
        //pion+
	if (pdgid > 0) {
	    charge = +1;
	} else {
	    charge = -1;
	}
	bib.Charge = charge;
	mass = 0.13957039;
	bib.Mass = mass;
    } else if (TMath::Abs(pdgid) == 2212) {
        //proton
	if (pdgid > 0) {
	    charge = +1;
	} else {
	    charge = -1;
	}
	bib.Charge = charge;
	mass = 0.93827208816;
	bib.Mass = mass;
    } else if (TMath::Abs(pdgid) == 2112) {
        //neutron
	charge = 0;
	bib.Charge = charge;
	mass = 0.93956542052;.
	bib.Mass = mass;
    } else {
	continue;
    }

    fOutputArray->Add(bib);
    fStableOutputArray->Add(bib);
  }
}
