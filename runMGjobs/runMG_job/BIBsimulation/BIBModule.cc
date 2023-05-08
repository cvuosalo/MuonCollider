/** \class BIBModule
 * 
 *  Psudo-simulation of Beam-Induced-Backgorund effect in multi-TeV Muon Collider with given distribution of momentum, position and particle type. Due to Delphes's modular structure, memory use of Delphes involving BIB simulation is proportional to the number of BIB particles required. DO NOT RUN WITH MORE THAN 10^5 PARICLES WITH LESS THAN 10GB RAM!!!
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

#include "TRandomGen.h"
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
    fThetaHistName = GetString("ThetaHistName", "theta");
    fPhiHistName = GetString("PhiHistName", "phi");
    fPdgIDHistName = GetString("PdgIDHistName", "pdgid");
    //fPdgEnergyHistName = GetString("PdgEnergyHistName", "pdgid_energy");
    fMomentumHistName = GetString("MomentumHistName", "px_py_pz");
 
    const char* mg5dir = std::getenv("mg5dir"); 
    fileName = GetString("FileName", "file");
    fileName = fileName.Insert(0, mg5dir);
    cout << "according to histograms from file " << fileName << "..." << endl;

    fInputArray = ImportArray(GetString("InputArray", "Delphes/allParticles"));
    fStableInputArray = ImportArray(GetString("StableInputArray", "Delphes/stableParticles"));
    //fItInputArray = fInputArray->MakeIterator();
    //fItStableInputArray = fStableInputArray->MakeIterator();

    fOutputArray = ExportArray(GetString("OutputArray", "allParticles"));
    fStableOutputArray = ExportArray(GetString("StableOutputArray", "stableParticles"));
    cout << endl << "Completing BIB Module initialization." << endl;
}

void BIBModule::Finish() {
    //if (fItInputArray) delete fItInputArray; 
    //if (fItStableInputArray) delete fItStableInputArray; 
}


void BIBModule::Process() {

    fItInputArray = fInputArray->MakeIterator();
    fItStableInputArray = fStableInputArray->MakeIterator();
    TFile* file = TFile::Open(fileName);
  
    TH1D* fxHist = (TH1D*) file -> Get(fxHistName);
    TH2D* fPositionHist = (TH2D*) file -> Get(fPositionHistName);
    TH3D* fMomentumHist = (TH3D*) file -> Get(fMomentumHistName);
    //TH2D* fPdgEnergyHist = (TH2D*) file -> Get(fPdgEnergyHistName);
    TH1D* fThetaHist = (TH1D*) file -> Get(fThetaHistName);
    TH1D* fPhiHist = (TH1D*) file -> Get(fPhiHistName);
    TH1D* fPdgIDHist = (TH1D*) file -> Get(fPdgIDHistName);
   /* 
    fxHist->SetBins(100,-600,600);
    fThetaHist->SetBins(100,-1.58,1.58);
    fPhiHist->SetBins(100,-1.58,1.58);
    fPdgIDHist->SetBins(100,0,2500);
*/
#if ROOT_VERSION_CODE < ROOT_VERSION(6,0,0)
    fxHist->SetBit(TH1::kCanRebin);
    fThetaHist->SetBit(TH1::kCanRebin);
    fPhiHist->SetBit(TH1::kCanRebin);
    fPdgIDHist->SetBit(TH1::kCanRebin);
#else
    fxHist->SetCanExtend(TH1::kXaxis);
    fThetaHist->SetCanExtend(TH1::kXaxis);
    fPhiHist->SetCanExtend(TH1::kXaxis);
    fPdgIDHist->SetCanExtend(TH1::kXaxis);
#endif
    Candidate *original;
    
    if (!fxHistName || !fPositionHistName || !fPhiHistName || !fThetaHistName || !fMomentumHistName || !fPdgIDHistName) return;

    fItInputArray -> Reset();
    fItStableInputArray -> Reset();

    while((original = static_cast<Candidate *>(fItInputArray -> Next()))) {
      original = static_cast<Candidate *>(original -> Clone());
      //cout << original -> Position.X()<< endl;
      fOutputArray -> Add(original);
    }


    while((original = static_cast<Candidate *>(fItStableInputArray -> Next()))) {
      original = static_cast<Candidate *>(original -> Clone());
      fStableOutputArray -> Add(original);
    }

    //cout << endl << "Completing original particle array copying." << endl;
    Candidate *bib;
    DelphesFactory *factory = GetFactory();
    Double_t pt, theta, eta, phi, px, py, pz, energy, mass, charge, pdgid, x, y, z, r;
    TLorentzVector bibPosition;
    TLorentzVector bibMomentum;
    	
    Int_t pdglist[7] = {11,13, 22, 111, 211, 2112, 2212};
    Int_t best = 0;


    gRandom->SetSeed(0);
    for (int i = 0; i < fNumParticles; i++) {
	    
        bib = factory -> NewCandidate();
	x = fxHist -> GetRandom();
    	fPositionHist -> GetRandom2(z, r);
	while (abs(x)>abs(r)) {
	    x = fxHist -> GetRandom();
	}
	fMomentumHist -> GetRandom3(px, py, pz);
	//fPdgEnergyHist -> GetRandom2(pdgid, energy);
	theta = fThetaHist -> GetRandom();
	while (abs(theta) > 1.570796) {
	    theta = fThetaHist -> GetRandom();
	}
	if (3.1415926-abs(theta)<=0.5) {
	    continue;
	}
	theta += 1.57079632;
	phi = (fPhiHist -> GetRandom())*2;
	
	pdgid = fPdgIDHist -> GetRandom();
	//cout <<  energy << endl ;
/* 
        if (energy <= 0) {
	    continue;
        }
  */
	eta = - TMath::Log(TMath::Tan(abs(theta) * 0.5));
	//if (abs(eta) > 10000) {
	  //  continue;
	//}
	//eta = (rnd.Rndm()) * 4 - 2;
	//pt = (rnd.Rndm()) * 10000; 
	pt = TMath::Sqrt(px*px + py*py);
	//cout << "pt is" << pt << endl;
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
	//cout << "x= " << x << " y= "<< y << " z= " << z<<endl;
	bibPosition.SetXYZT(x, y, z, 0);
	//bibMomentum.SetPxPyPzE(px, py, pz, energy);
    
	bib -> PID = pdgid; 
	bib -> D1 = -1;
	bib -> D2 = -1;
	bib -> M1 = -1;
	bib -> M2 = -1;
	bib -> Status = 1; 
	bib -> IsPU = 0;

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
        //energy = TMath::Sqrt(px*px+py*py+pz*pz+mass*mass); 
	bibMomentum.SetPtEtaPhiM(pt, eta, phi, mass);
	//bibMomentum.SetPxPyPzE(px, py, pz, energy);

	bib -> Position = bibPosition;
	bib -> Momentum = bibMomentum;
	
	fOutputArray -> Add(bib);
	fStableOutputArray -> Add(bib);

        if (bibMomentum.E() <= 0) {
	    continue;
        }
	
	//cout << endl << "Not crashing in " << i <<"th BIB module loops." << endl;
    }
   
    //cout << endl << "counting" << allcounter << "all particles";
    //cout << endl << "counting" << stablecounter << "stable particles" ;
    //cout << endl << "counting" << BIBcounter << "BIB particles added" << endl;

    delete original;
    delete fxHist;
    delete fPositionHist;
    delete fPdgIDHist;
    delete fThetaHist;
    delete fPhiHist;
    delete fMomentumHist;

    file -> Close();
    if (file) delete file;
    //cout << endl << "BIB particles successfully adding to event." << endl;

    if (fItInputArray) delete fItInputArray; 
    if (fItStableInputArray) delete fItStableInputArray; 

}
