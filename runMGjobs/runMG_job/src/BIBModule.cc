/** \class BIBModule
 
 * 
 
 *  Psudo-simulation of BIB effect with given distribution of energy and position.
 
 * 
 
 *  \author H. Jia - UW-Madison
 
 *
 
 */

#include "modules/BIBModule.h"

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
#include "TMath.h"

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
/*
  TFile* file = TFile::Open(fileName);
  
  TH1D* fxHist = (TH1D*) file->Get(fxHistName);
  TH2D* fPositionHist = (TH2D*) file->Get(fPositionHistName);
  TH1D* fPhiHist = (TH1D*) file->Get(fPhiHistName);
  TH1D* fThetaHist = (TH1D*) file->Get(fThetaHistName);
  TH2D* fPdgEnergyHist = (TH2D*) file->Get(fPdgEnergyHistName);
  TH3D* fMomentumHist = (TH3D*) file->Get(fMomentumHistName);
*/
}

void BIBModule::Finish()
{
  if(fItInputArray) delete fItInputArray; 
  if(fItStableInputArray) delete fItStableInputArray; 
}

void BIBModule::Process()
{

  TFile* file = TFile::Open(fileName);
  
  fxHist = (TH1D*) file->Get(fxHistName);
  fPositionHist = (TH2D*) file->Get(fPositionHistName);
  fPhiHist = (TH1D*) file->Get(fPhiHistName);
  fThetaHist = (TH1D*) file->Get(fThetaHistName);
  fPdgEnergyHist = (TH2D*) file->Get(fPdgEnergyHistName);
  fMomentumHist = (TH3D*) file->Get(fMomentumHistName);

  GenParticle *original;
  GenParticle *bib = new GenParticle;
  Double_t px, py, pz, energy, theta, eta, phi, mass, charge, pdgid, x, y, z, r;

  if (!fxHistName) return;
  if (!fPositionHistName) return;
  if (!fPhiHistName) return;
  if (!fThetaHistName) return;
  if (!fMomentumHistName) return;
  if (!fPdgEnergyHistName) return;

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

  for (Int_t i = 0; i < fNumParticles; i++)
  {
    phi = fPhiHist -> GetRandom();
    theta = fThetaHist -> GetRandom();
    x = fxHist -> GetRandom();
    fPositionHist -> GetRandom2(z, r);
    fMomentumHist -> GetRandom3(px, py, pz);
    fPdgEnergyHist -> GetRandom2(pdgid, energy);
    
    if(energy <= 0.0) continue;
    
    bib -> X = x;
    y = TMath::Sqrt(r * r - x * x);
    bib -> Y = y;
    bib -> Z = z;
    eta = - TMath::Log(TMath::Tan(theta/2));
    bib -> Eta = eta;
    bib -> Phi = phi;
    bib -> E = energy;
    bib -> Px = px;
    bib -> Py = py;
    bib -> Pz = pz;

    bib -> PID = pdgid; 

    //assign mass and charge by pdgID here

    if (TMath::Abs(pdgid) == 11) {
        //electron
	if (pdgid > 0) {
	    charge = -1;
	} else {
	    charge = +1;
	}
	bib -> Charge = charge;
	mass = 0.00051099895;
	bib -> Mass = mass;
    } else if (TMath::Abs(pdgid) == 13) {
        //muon
	if (pdgid > 0) {
	    charge = -1;
	} else {
	    charge = +1;
	}
	charge = -1;
	bib -> Charge = charge;
	mass = 0.1056583755;
	bib -> Mass = mass;
    } else if (TMath::Abs(pdgid) == 22) {
	//photon
	charge = 0;
	bib -> Charge = charge;
	mass = 0;
	bib -> Mass = mass;
    } else if (TMath::Abs(pdgid) == 111) {
        //neutral pion
	charge = 0;
	bib -> Charge = charge;
	mass = 0.1349768;
	bib -> Mass = mass;
    } else if (TMath::Abs(pdgid) == 211) {
        //pion+
	if (pdgid > 0) {
	    charge = +1;
	} else {
	    charge = -1;
	}
	bib -> Charge = charge;
	mass = 0.13957039;
	bib -> Mass = mass;
    } else if (TMath::Abs(pdgid) == 2212) {
        //proton
	if (pdgid > 0) {
	    charge = +1;
	} else {
	    charge = -1;
	}
	bib -> Charge = charge;
	mass = 0.93827208816;
	bib -> Mass = mass;
    } else if (TMath::Abs(pdgid) == 2112) {
        //neutron
	charge = 0;
	bib -> Charge = charge;
	mass = 0.93956542052;
	bib -> Mass = mass;
    } else {
	continue;
    }

    fOutputArray->Add(bib);
    fStableOutputArray->Add(bib);
  }
}
