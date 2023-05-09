/** \class BIBECal
 * 
 *  Psudo-simulation of BIB effect with given distribution of energy and position on the electromagnetic calorimeter. Only charged hadrons, photons, and leptons in proximity of original EnergyFlow photons and EnergyFlow tracks are considered. Particle propagation of BIB is built-in in this module.
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
#include "TFile.h"
#include "TTree.h"

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
    //Aquire input parameters
    fNumParticles = GetInt("NumParticles", 1000000);
    fPhotonsDeltaR = GetDouble("PhotonsDeltaR", 0.5);
    fTracksDeltaR = GetDouble("TracksDeltaR", 0.5);

    fBz = GetDouble("Bz", 0.0);

    //cout << endl << "Generating " << fNumParticles << " BIB particles ";
    fxHistName = GetString("xHistName", "x");
    fPositionHistName = GetString("PositionHistName", "z_r");
    fThetaHistName = GetString("ThetaHistName", "theta");
    fPhiHistName = GetString("PhiHistName", "phi");
    fPdgIDHistName = GetString("PdgIDHistName", "pdgid");
    fMomentumHistName = GetString("MomentumHistName", "px_py_pz");
 
    const char* mg5dir = std::getenv("mg5dir"); 
    fileName = GetString("FileName", "file");
    fileName = fileName.Insert(0, mg5dir);
    //cout << "according to histograms from file " << fileName << "..." << endl;

    fPhotonsInputArray = ImportArray(GetString("PhotonsInputArray", ""));
    fTracksInputArray = ImportArray(GetString("TracksInputArray", ""));
    fItPhotonsInputArray = fPhotonsInputArray->MakeIterator();
    fItTracksInputArray = fTracksInputArray->MakeIterator();

    fPhotonsOutputArray = ExportArray(GetString("PhotonsOutputArray", "eflowPhotons"));
    fTracksOutputArray = ExportArray(GetString("TracksOutputArray", "eflowTracks"));

    fItPhotonsOutputArray = fPhotonsOutputArray->MakeIterator();
    fItTracksOutputArray = fTracksOutputArray->MakeIterator();

    //cout << endl << "Completing BIB Module initialization." << endl;
}

void BIBECal::Finish() {
    //Clean up memory use
    if (fItPhotonsInputArray) delete fItPhotonsInputArray; 
    if (fItTracksInputArray) delete fItTracksInputArray; 
    if (fItPhotonsOutputArray) delete fItPhotonsOutputArray; 
    if (fItTracksOutputArray) delete fItTracksOutputArray; 
}


void BIBECal::Process() {
    //Get Full Simulation distribution from given root file
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

    fItPhotonsInputArray -> Reset();
    fItTracksInputArray -> Reset();

    while((original = static_cast<Candidate *>(fItPhotonsInputArray -> Next()))) {
      original = static_cast<Candidate *>(original -> Clone());
      fPhotonsOutputArray -> Add(original);
    }


    while((original = static_cast<Candidate *>(fItTracksInputArray -> Next()))) {
      original = static_cast<Candidate *>(original -> Clone());
      fTracksOutputArray -> Add(original);
    }

    // Generate randomize particles according the distribution
    Double_t pt, theta, eta, phi, px, py, pz, mass, charge, pdgid, x, y, z, r;
    TLorentzVector bibMomentum;
    Double_t gammam, e, omega, r_helix;
    Double_t x_c, y_c, td, pio, phid, phi_0;
    const Double_t c_light = 2.99792458E8;
    	
    Int_t pdglist[7] = {11,13, 22, 111, 211, 2112, 2212};
    Int_t best = 0;

    bool EMFlag;
    bool HadFlag;
    bool propaFlag;

    gRandom->SetSeed(0);

    for (int i = 0; i < fNumParticles; ++i) {
	EMFlag = false;
	HadFlag = false;
	propaFlag = false;
	    
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

	//Check if PdgID is correct and gives corresponding charge and mass
	switch((Int_t) TMath::Abs(pdgid)) {
	    case 11:
		mass = 0.00051099895;
		charge = (pdgid > 0) ? -1 : 1;
		EMFlag = true;
		propaFlag = true;
		break;
	    case 13:
		mass = 0.1056583755;
		charge = (pdgid > 0) ? -1 : 1;
		EMFlag = true;
		propaFlag = true;
		break;
	    case 22:
		mass = 0;
		charge = 0;
		EMFlag = true;
		break;
	    case 111:
		continue;
		mass = 0.1349768;
		charge = 0;
		break;
	    case 211:
		HadFlag = true;
		mass = 0.13957039;
		charge = (pdgid > 0) ? -1 : 1;
		propaFlag = true;
		break;
	    case 2212:
		HadFlag = true;
		mass = 0.93827208816;
		charge = (pdgid > 0) ? -1 : 1;
		propaFlag = true;
		break;
	    case 2112:
		continue;
		mass = 0.93956542052;
		charge = 0;
		break;
	    default:
		cout << "Error: Not a valid PDG ID!" << endl;
		continue;
		break;
	}
	
	x = fxHist -> GetRandom();
    	fPositionHist -> GetRandom2(z, r);
	while (abs(x)>abs(r)) {
	    x = fxHist -> GetRandom();
	}
	y = TMath::Sqrt(r * r - x * x);
	
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

	//Particle Propagation for Charge Particles

	if (propaFlag) {
	    // initial transverse momentum p_{T0}: Part->pt
	    // initial transverse momentum direction phi_0 = -atan(p_{X0} / p_{Y0})
	    // relativistic gamma: gamma = E / mc^2; gammam = gamma * m
	    // gyration frequency omega = q * Bz / (gammam)
	    // helix radius r = p_{T0} / (omega * gammam)
	     
	    e = bibMomentum.E();
	    px = bibMomentum.Px();
	    py = bibMomentum.Py();
	    gammam = e * 1.0E9 / (c_light * c_light); // gammam in [eV/c^2] 
	    omega = charge * fBz / gammam; // omega is here in [89875518/s]
	    r_helix = pt / (charge * fBz) * 1.0E9 / c_light; // in [m]

	    phi_0 = TMath::ATan2(py, px); // [rad] in [-pi, pi]

	    //helix axis coordinates
	    x_c = x + r_helix * TMath::Sin(phi_0);
	    y_c = y - r_helix * TMath::Cos(phi_0);
	    
	    // time of closest approach
	    td = (phi_0 + TMath::ATan2(x_c, y_c)) / omega;

	    // remove all the modulo pi that might have come from the atan
	    pio = TMath::Abs(TMath::Pi() / omega);
	    while(TMath::Abs(td) > 0.5 * pio) {
		td -= TMath::Sign(1.0, td) * pio;
	    }
	    phid = phi_0 - omega * td;

	    // momentum at closest approach
	    bibMomentum.SetPtEtaPhiM(pt, eta, phid, mass);
	}

	if (fItPhotonsOutputArray) delete fItPhotonsOutputArray; 
        if (fItTracksOutputArray) delete fItTracksOutputArray; 

        fItPhotonsOutputArray = fPhotonsOutputArray->MakeIterator();
        fItTracksOutputArray = fTracksOutputArray->MakeIterator();

	fItPhotonsInputArray -> Reset();
	fItTracksInputArray -> Reset();

	//check if BIB generated in proximity of original eflowPhotons
	//bool EMaddornot = false;
	//bool HADaddornot = false;
	if (EMFlag) {
	    TLorentzVector originalMomentum;
	    while((original = static_cast<Candidate *>(fItPhotonsOutputArray -> Next()))) {
		originalMomentum = original -> Momentum;
		if (originalMomentum.DeltaR(bibMomentum) < fPhotonsDeltaR){
		    originalMomentum = originalMomentum + bibMomentum;
		    original -> Momentum = originalMomentum;
		} 
            }
	}
	

	if (HadFlag) {
	    TLorentzVector originalMomentum;
	    while((original = static_cast<Candidate *>(fItTracksOutputArray -> Next()))) {
		originalMomentum = original -> Momentum;
		if (originalMomentum.DeltaR(bibMomentum) < fTracksDeltaR){
		    originalMomentum = originalMomentum + bibMomentum;
		    original -> Momentum = originalMomentum;
		} 
            }
	}
	

	if (fItPhotonsInputArray) delete fItPhotonsInputArray; 
        if (fItTracksInputArray) delete fItTracksInputArray; 

        fItPhotonsInputArray = fPhotonsInputArray->MakeIterator();
        fItTracksInputArray = fTracksInputArray->MakeIterator();

	fItPhotonsInputArray -> Reset();
	fItTracksInputArray -> Reset();

	//Calculate ratio of pT after adding BIB 
	Int_t ph=0;
	Int_t tr=0;
	Double_t Photons[fPhotonsInputArray -> GetEntries()]; 
	Double_t Tracks[fTracksInputArray -> GetEntries()]; 
	while((original = static_cast<Candidate *>(fItPhotonsInputArray -> Next()))) {
	    TLorentzVector originalMomentum;
	    originalMomentum = original -> Momentum;
	    Photons[ph] = 1/originalMomentum.Pt();
	    ph++;
	}

	while((original = static_cast<Candidate *>(fItTracksInputArray -> Next()))) {
	    TLorentzVector originalMomentum;
	    originalMomentum = original -> Momentum;
	    Tracks[tr] = 1/originalMomentum.Pt();
	    tr++;
	}
	
	if (fItPhotonsOutputArray) delete fItPhotonsOutputArray; 
        if (fItTracksOutputArray) delete fItTracksOutputArray; 

        fItPhotonsOutputArray = fPhotonsOutputArray->MakeIterator();
        fItTracksOutputArray = fTracksOutputArray->MakeIterator();

	fItPhotonsOutputArray -> Reset();
	fItTracksOutputArray -> Reset();

	ph = 0;
	tr = 0;
	while((original = static_cast<Candidate *>(fItPhotonsOutputArray -> Next()))) {
	    TLorentzVector originalMomentum;
	    originalMomentum = original -> Momentum;
	    original -> dNdx = Photons[ph] * originalMomentum.Pt();
	    ph++;
	}
	while((original = static_cast<Candidate *>(fItTracksOutputArray -> Next()))) {
	    TLorentzVector originalMomentum;
	    originalMomentum = original -> Momentum;
	    original -> dNdx = Tracks[tr] * originalMomentum.Pt();
	    tr++;
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
