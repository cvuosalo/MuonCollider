/** \class BIBModule
 * 
 *  Psudo-simulation of BIB effect with given distribution of energy and position.
 * 
 *  \author H. Jia - UW-Madison
 *
 */

#include "modules/BIBNeutralHadrons.h"

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

#include "TRandomGen.h"
#include "TMath.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

//------------------------------------------------------------------------------

BIBNeutralHadrons::BIBNeutralHadrons() : 
    fInputArray(0)
{
}

//------------------------------------------------------------------------------

BIBNeutralHadrons::~BIBNeutralHadrons() 
{
}

//------------------------------------------------------------------------------

void BIBNeutralHadrons::Init()
{
    //Get parameters
    fNumParticles = GetInt("NumParticles", 1000000);
    fDeltaR = GetDouble("DeltaR", 0.1);
    cout << endl << "Generating " << fNumParticles << " BIB particles ";

    fxHistName = GetString("xHistName", "x");
    fPositionHistName = GetString("PositionHistName", "z_r");
    fThetaHistName = GetString("ThetaHistName", "theta");
    fPhiHistName = GetString("PhiHistName", "phi");
    fPdgIDHistName = GetString("PdgIDHistName", "pdgid");
    fMomentumHistName = GetString("MomentumHistName", "px_py_pz");
 
    const char* mg5dir = std::getenv("mg5dir"); 
    fileName = GetString("FileName", "file");
    fileName = fileName.Insert(0, mg5dir);
    cout << "according to histograms from file " << fileName << "..." << endl;

    fInputArray = ImportArray(GetString("InputArray", ""));
    fItInputArray = fInputArray->MakeIterator();

    fOutputArray = ExportArray(GetString("OutputArray", ""));
    fItOutputArray = fOutputArray->MakeIterator();
    cout << endl << "Completing BIB Module initialization." << endl;
}

void BIBNeutralHadrons::Finish() {
    //Clean up memory use
    if (fItInputArray) delete fItInputArray; 
    if (fItOutputArray) delete fItOutputArray; 
}


void BIBNeutralHadrons::Process() {
    //Acquire Full simulation distribution from given root file
    TFile* file = TFile::Open(fileName);
  
    TH1D* fxHist = (TH1D*) file -> Get(fxHistName);
    TH2D* fPositionHist = (TH2D*) file -> Get(fPositionHistName);
    TH3D* fMomentumHist = (TH3D*) file -> Get(fMomentumHistName);
    TH1D* fThetaHist = (TH1D*) file -> Get(fThetaHistName);
    TH1D* fPhiHist = (TH1D*) file -> Get(fPhiHistName);
    TH1D* fPdgIDHist = (TH1D*) file -> Get(fPdgIDHistName);
   
    //Copy original array
    Candidate *original;
    
    if (!fxHistName || !fPositionHistName || !fPhiHistName || !fThetaHistName || !fMomentumHistName || !fPdgIDHistName) return;

    fItInputArray -> Reset();

    while((original = static_cast<Candidate *>(fItInputArray -> Next()))) {
      original = static_cast<Candidate *>(original -> Clone());
      fOutputArray -> Add(original);
    }

    Double_t pt, theta, eta, phi, px, py, pz, energy, mass, charge, pdgid, x, y, z, r;
    TLorentzVector bibMomentum;
    	
    Int_t pdglist[7] = {11,13, 22, 111, 211, 2112, 2212};
    Int_t best = 0;

    bool proximity;
    bool HadFlag;

    gRandom->SetSeed(0);

    for (int i = 0; i < fNumParticles; ++i) {
	HadFlag = false;
	proximity = false;
	    
	pdgid = fPdgIDHist -> GetRandom();

	for (int i=1; i < 7; i++) {
      	    if (abs(abs(pdgid)-pdglist[i]) < abs(abs(pdgid)-pdglist[best])) {
        	best = i;
	    } 
	}

        if (pdgid > 0) {
	    pdgid = pdglist[best];
	} else {
            pdgid = -pdglist[best];
	}

	//Check PdgID and give corresponding charge and mass
	switch((Int_t) TMath::Abs(pdgid)) {
	    case 11:
		continue;
		mass = 0.00051099895;
		charge = (pdgid > 0) ? -1 : 1;
		break;
	    case 13:
		continue;
		mass = 0.1056583755;
		charge = (pdgid > 0) ? -1 : 1;
		break;
	    case 22:
		continue;
		mass = 0;
		charge = 0;
		break;
	    case 111:
		HadFlag = true;
		mass = 0.1349768;
		charge = 0;
		break;
	    case 211:
		continue;
		mass = 0.13957039;
		charge = (pdgid > 0) ? -1 : 1;
		break;
	    case 2212:
		continue;
		mass = 0.93827208816;
		charge = (pdgid > 0) ? -1 : 1;
		break;
	    case 2112:
		HadFlag = true;
		mass = 0.93956542052;
		charge = 0;
		break;
	    default:
		cout << "Error: Not a valid PDG ID!" << endl;
		continue;
		break;
	}
        
	fMomentumHist -> GetRandom3(px, py, pz);
	theta = fThetaHist -> GetRandom();
	while (abs(theta) > 1.570796) {
	    theta = fThetaHist -> GetRandom();
	}
	if (3.1415926-abs(theta)<=0.5) {
	    continue;
	}
	theta += 1.57079632;
	phi = (fPhiHist -> GetRandom())*2;
	
	eta = - TMath::Log(TMath::Tan(abs(theta) * 0.5));

	pt = TMath::Sqrt(px*px + py*py);
	

	bibMomentum.SetPtEtaPhiM(pt, eta, phi, mass);

	if (fItOutputArray) delete fItOutputArray; 
	fItOutputArray = fOutputArray -> MakeIterator();
	fItOutputArray -> Reset();

	//check if BIB generated in proximity of original eflowPhotons

	if (HadFlag) {
	    TLorentzVector originalMomentum;
	    while((original = static_cast<Candidate *>(fItOutputArray -> Next()))) {
		originalMomentum = original -> Momentum;
		if (originalMomentum.DeltaR(bibMomentum) < fDeltaR){
		    originalMomentum = originalMomentum + bibMomentum;
		    original -> Momentum = originalMomentum;
		} 
            }
	}

        if (fItInputArray) delete fItInputArray; 

        fItInputArray = fInputArray->MakeIterator();

	fItInputArray -> Reset();

	Int_t had=0;
	Double_t neuHads[fInputArray -> GetEntries()]; 
	while((original = static_cast<Candidate *>(fItInputArray -> Next()))) {
	    TLorentzVector originalMomentum;
	    originalMomentum = original -> Momentum;
	    neuHads[had] = 1/originalMomentum.Pt();
	    had++;
	}
	
        if (fItOutputArray) delete fItOutputArray; 

        fItOutputArray = fOutputArray->MakeIterator();

	fItOutputArray -> Reset();

	had = 0;
	
	while((original = static_cast<Candidate *>(fItOutputArray -> Next()))) {
	    TLorentzVector originalMomentum;
	    originalMomentum = original -> Momentum;
	    original -> dNdx = neuHads[had] * originalMomentum.Pt();
	    had++;
	}
    }

    delete original;
    delete fxHist;
    delete fPositionHist;
    delete fPdgIDHist;
    delete fThetaHist;
    delete fPhiHist;
    delete fMomentumHist;

    file -> Close();
    if (file) delete file;
    
}
