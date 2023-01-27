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

    cout << endl << "Generating " << fNumParticles << " BIB particles ";
    fxHistName = GetString("xHistName", "x");
    fPositionHistName = GetString("PositionHistName", "z_r");
    //fPhiHistName = GetString("PhiHistName", "phi");
    //fThetaHistName = GetString("ThetaHistName", "theta");
    fPdgEnergyHistName = GetString("PdgEnergyHistName", "pdgid_energy");
    fMomentumHistName = GetString("MomentumHistName", "px_py_pz");
  
    fileName = GetString("FileName", "file");
    cout << "according to histograms from file " << fileName << "..." << endl;

    fInputArray = ImportArray(GetString("InputArray", "Delphes/allParticles"));
    fStableInputArray = ImportArray(GetString("StableInputArray", "Delphes/stableParticles"));
    fItInputArray = fInputArray->MakeIterator();
    fItStableInputArray = fStableInputArray->MakeIterator();

    fOutputArray = ExportArray(GetString("OutputArray", "allParticles"));
    fStableOutputArray = ExportArray(GetString("StableOutputArray", "stableParticles"));
 /*  
    file = TFile::Open(fileName);
  
    fxHist = (TH1D*) file -> Get(fxHistName);
    fPositionHist = (TH2D*) file -> Get(fPositionHistName);
    //fPhiHist = (TH1D*) file -> Get(fPhiHistName);
    //fThetaHist = (TH1D*) file -> Get(fThetaHistName);
    fPdgEnergyHist = (TH2D*) file -> Get(fPdgEnergyHistName);
    fMomentumHist = (TH3D*) file -> Get(fMomentumHistName);
*/
    //Candidate *bib = new Candidate;
}

void BIBModule::Finish() {
    if (fItInputArray) delete fItInputArray; 
    if (fItStableInputArray) delete fItStableInputArray; 
}

void BIBModule::Process() {

    TFile* file = TFile::Open(fileName);
  
    TH1D* fxHist = (TH1D*) file -> Get(fxHistName);
    TH2D* fPositionHist = (TH2D*) file -> Get(fPositionHistName);
    //TH1D* fPhiHist = (TH1D*) file -> Get(fPhiHistName);
    //TH1D* fThetaHist = (TH1D*) file -> Get(fThetaHistName);
    TH2D* fPdgEnergyHist = (TH2D*) file -> Get(fPdgEnergyHistName);
    TH3D* fMomentumHist = (TH3D*) file -> Get(fMomentumHistName);
    

    Candidate *original;
    //std::unique_ptr<Candidate> bib(new Candidate());
    
    //TLorentzVector bibPosition, bibMomentum;
    //Double_t  px, py, pz, energy, mass, charge, pdgid, x, y, z, r;
/*
    Int_t allcounter, stablecounter, BIBcounter;*/
    //allcounter = 0;
    //stablecounter = 0;
    //BIBcounter = 0;

    if (!fxHistName || !fPositionHistName || /*!fPhiHistName || !fThetaHistName ||*/ !fMomentumHistName || !fPdgEnergyHistName) return;

    fItInputArray -> Reset();
    fItStableInputArray -> Reset();

    while((original = static_cast<Candidate *>(fItInputArray -> Next()))) {
      original = static_cast<Candidate *>(original -> Clone());
      fOutputArray -> Add(original);
      //allcounter++;
    }


    while((original = static_cast<Candidate *>(fItStableInputArray -> Next()))) {
      original = static_cast<Candidate *>(original -> Clone());
      fStableOutputArray -> Add(original);
      //stablecounter++;
    }

    
    
    
    for (int i = 0; i < fNumParticles; i++) {
	Candidate *bib;
	DelphesFactory *factory = GetFactory();
	Double_t px, py, pz, energy, mass, charge, pdgid, x, y, z, r;
	TLorentzVector bibPosition;
	TLorentzVector bibMomentum;
    
	bib = factory -> NewCandidate();
	//phi = fPhiHist -> GetRandom();
	//theta = fThetaHist -> GetRandom();
	x = fxHist -> GetRandom();
    	fPositionHist -> GetRandom2(z, r);
	fMomentumHist -> GetRandom3(px, py, pz);
	fPdgEnergyHist -> GetRandom2(pdgid, energy);
 
        if (energy <= 0) {
	    continue;
        }
  
	//eta = - TMath::Log(TMath::Tan(theta/2));
	//pt = TMath::Sqrt(px*px + py*py);
	
	Int_t pdglist[7] = {11,13, 22, 111, 211, 2112, 2212};
	Int_t best = 0;

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

        y = TMath::Sqrt(r * r - x * x);
	bibPosition.SetXYZT(x, y, z, 0);
	bibMomentum.SetPxPyPzE(px, py, pz, energy);
	//bibMomentum.SetPtEtaPhiE(pt, eta, phi, energy);
    
	bib -> PID = pdgid; 
	bib -> D1 = -1;
	bib -> D2 = -1;
	bib -> M1 = -1;
	bib -> M2 = -1;
	bib -> Status = 1; 
	bib -> IsPU = 0;

	bib -> Position = bibPosition;
	bib -> Momentum = bibMomentum;
	
	//assign mass and charge by pdgID here
	switch((Int_t) TMath::Abs(pdgid)) {
	    case 11:
		mass = 0.00051099895;
		charge = (pdgid > 0) ? -1 : 1;
		break;
	    case 13:
		mass = 0.1056583755;
		charge = (pdgid > 0) ? -1 : 1;
		break;
	    case 22:
		mass = 0;
		charge = 0;
		break;
	    case 111:
		mass = 0.1349768;
		charge = 0;
		break;
	    case 211:
		mass = 0.13957039;
		charge = (pdgid > 0) ? -1 : 1;
		break;
	    case 2212:
		mass = 0.93827208816;
		charge = (pdgid > 0) ? -1 : 1;
		break;
	    case 2112:
		mass = 0.93956542052;
		charge = 0;
		break;
	    default:
		cout << "Error: Not a valid PDG ID!" << endl;
		continue;
		break;
	} 
 
	bib -> Charge = charge;
	bib -> Mass = mass;
    
	fOutputArray -> Add(bib);
	fStableOutputArray -> Add(bib);

        if (bibMomentum.E() <= 0) {
	    continue;
        }
    }
   
    //cout << endl << "counting" << allcounter << "all particles";
    //cout << endl << "counting" << stablecounter << "stable particles" ;
    //cout << endl << "counting" << BIBcounter << "BIB particles added" << endl;

    delete original;
    //if(bib) delete bib;
    delete fxHist;
    delete fPositionHist;
    //if(fPhiHist) delete fPhiHist;
    //if(fThetaHist) delete fThetaHist;
    delete fPdgEnergyHist;

    file -> Close();
    if (file) delete file;
}
