/** \class BIBModule
 * 
 *  Psudo-simulation of BIB effect with given distribution of energy and position.
 * 
 *  \author H. Jia - UW-Madison
 *
 */

#include "modules/BIBECal.h"

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

BIBECal::BIBECal() : 
    fPhotonsInputArray(0), fTracksInputArray(0)
{
}

//------------------------------------------------------------------------------

BIBECal::~BIBECal() 
{
}

//------------------------------------------------------------------------------

void BIBECal::Init()
{
    fNumParticles = GetInt("NumParticles", 30000000);
    fPhotonsDeltaR = GetDouble("PhotonsDeltaR", 0.5);
    fTracksDeltaR = GetDouble("TracksDeltaR", 0.5);

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

    fPhotonsInputArray = ImportArray(GetString("PhotonsInputArray", ""));
    fStableInputArray = ImportArray(GetString("TracksInputArray", ""));
    //fItInputArray = fInputArray->MakeIterator();
    //fItStableInputArray = fStableInputArray->MakeIterator();

    fPhotonsOutputArray = ExportArray(GetString("PhotonsOutputArray", "eflowPhotons"));
    fTracksOutputArray = ExportArray(GetString("TracksOutputArray", "eflowTracks"));
    cout << endl << "Completing BIB Module initialization." << endl;
}

void BIBECal::Finish() {
    //if (fItInputArray) delete fItInputArray; 
    //if (fItStableInputArray) delete fItStableInputArray; 
}


void BIBECal::Process() {

    fItPhotonsInputArray = fPhotonsInputArray->MakeIterator();
    fItTracksInputArray = fTracksInputArray->MakeIterator();
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

    fItPhotonsInputArray -> Reset();
    fItTracksInputArray -> Reset();

    TObjArray *originalPhotons;
    TObjArray *originalTracks;
    while((original = static_cast<Candidate *>(fItPhotonsInputArray -> Next()))) {
      original = static_cast<Candidate *>(original -> Clone());
      //cout << original -> Position.X()<< endl;
      fPhotonsOutputArray -> Add(original);
      originalPhotons -> Add(original);
    }


    while((original = static_cast<Candidate *>(fItTracksInputArray -> Next()))) {
      original = static_cast<Candidate *>(original -> Clone());
      fTracksOutputArray -> Add(original);
      originalTracks -> Add(original);
    }

    //cout << endl << "Completing original particle array copying." << endl;
    Candidate *bib;
    DelphesFactory *factory = GetFactory();
    Double_t pt, theta, eta, phi, px, py, pz, energy, mass, charge, pdgid, x, y, z, r;
    TLorentzVector bibPosition;
    TLorentzVector bibMomentum;
    	
    Int_t pdglist[7] = {11,13, 22, 111, 211, 2112, 2212};
    Int_t best = 0;

    bool EMFlag;
    bool proximity;
    bool HadFlag;

    gRandom->SetSeed(0);
    for (int i = 0; i < fNumParticles; i++) {
	EMFlag = false;
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

	switch((Int_t) TMath::Abs(pdgid)) {
	    case 11:
		mass = 0.00051099895;
		charge = (pdgid > 0) ? -1 : 1;
		EMFlag = true;
		break;
	    case 13:
		mass = 0.1056583755;
		charge = (pdgid > 0) ? -1 : 1;
		EMFlag = true;
		break;
	    case 22:
		mass = 0;
		charge = 0;
		EMFlag = true;
		break;
	    case 111:
		continue;
		//HadFlag = true;
		mass = 0.1349768;
		charge = 0;
		break;
	    case 211:
		HadFlag = true;
		mass = 0.13957039;
		charge = (pdgid > 0) ? -1 : 1;
		break;
	    case 2212:
		HadFlag = true;
		mass = 0.93827208816;
		charge = (pdgid > 0) ? -1 : 1;
		break;
	    case 2112:
		continue;
		//HadFlag = true;
		mass = 0.93956542052;
		charge = 0;
		break;
	    default:
		cout << "Error: Not a valid PDG ID!" << endl;
		continue;
		break;
	}
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
 
	bib -> Charge = charge;
	bib -> Mass = mass;
        //energy = TMath::Sqrt(px*px+py*py+pz*pz+mass*mass); 
	bibMomentum.SetPtEtaPhiM(pt, eta, phi, mass);
	//bibMomentum.SetPxPyPzE(px, py, pz, energy);

	bib -> Position = bibPosition;
	bib -> Momentum = bibMomentum;

	//check if BIB generated in proximity of original eflowPhotons
	if (EMFlag) {
	    while((original = static_cast<Candidate *>(originalPhotons -> Next()))) {
     	        original = static_cast<Candidate *>(original -> Clone());
		TLorentzVector *originalMomentum = original -> Momentum;
		if (originalMomentum.DeltaR(bibMomentum) < fPhotonsDeltaR) proximity = true;
            }
	    if (proximity) fPhotonsOutputArray -> Add(bib);
	}

	if (HadFlag) {
	    while((original = static_cast<Candidate *>(originalTracks -> Next()))) {
     	        original = static_cast<Candidate *>(original -> Clone());
		TLorentzVector *originalMomentum = original -> Momentum;
		if (originalMomentum.DeltaR(bibMomentum) < fTracksDeltaR) proximity = true;
            }
	    if (proximity) fTracksOutputArray -> Add(bib);
	}
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

    if (fItPhotonsInputArray) delete fItPhotonsInputArray; 
    if (fItTracksInputArray) delete fItTracksInputArray; 

}
