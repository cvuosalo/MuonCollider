#ifdef __CLING__

R__LOAD_LIBRARY(libDelphes)

#include "classes/DelphesClasses.h"

#include "classes/DelphesFactory.h"

#include "classes/DelphesStream.h"

#include "classes/SortableObject.h"

#include "modules/Delphes.h"

#include "external/ExRootAnalysis/ExRootProgressBar.h"

#include "external/ExRootAnalysis/ExRootTreeBranch.h"

#include "external/ExRootAnalysis/ExRootTreeReader.h"

#include "external/ExRootAnalysis/ExRootTreeWriter.h"

#include "external/ExRootAnalysis/ExRootTask.h"

#include <string.h>

#endif

#include <TVector.h>

#include <TTree.h>

#include <TBranch.h>

#include <TLeaf.h>

#include <TH1.h>

#include <TH2.h>

#include <TCanvas.h>

#include <TMath.h>

#include <TLorentzVector.h>

#include <TLatex.h>

#include <time.h>

void JetEnergyFix(Double_t AKTjeteta,
    Double_t AKTjetPt,
	Double_t JER[10][10]) { // Double_t MuonwJER[10][10], Double_t MuonwoJER[10][10], bool MuonTagging){
    
    Double_t AKTjetTheta = 2 *atan(exp(-AKTjeteta));
    Int_t ThetaGrid;
    Int_t PTGrid;

    if (0 <= AKTjetTheta and AKTjetTheta <= 0.3) {
        ThetaGrid = 0;
    }
    if (0.3 < AKTjetTheta and AKTjetTheta <= 0.6) {
        ThetaGrid = 1;
    }
    if (0.6 < AKTjetTheta and AKTjetTheta <= 0.9) {
        ThetaGrid = 2;
    }
    if (0.9 < AKTjetTheta and AKTjetTheta <= 1.2) {
        ThetaGrid = 3;
    }
    if (1.2 < AKTjetTheta and AKTjetTheta <= 1.5) {
        ThetaGrid = 4;
    }
    if (1.5 < AKTjetTheta and AKTjetTheta <= 1.8) {
        ThetaGrid = 5;
    }
    if (1.8 < AKTjetTheta and AKTjetTheta <= 2.1) {
        ThetaGrid = 6;
    }
    if (2.1 < AKTjetTheta and AKTjetTheta <= 2.4) {
        ThetaGrid = 7;
    }
    if (2.4 < AKTjetTheta and AKTjetTheta <= 2.7) {
        ThetaGrid = 8;
    }
    if (2.7 < AKTjetTheta) {
        ThetaGrid = 9;
    }
    if (0 < AKTjetPt and AKTjetPt <= 50) {
        PTGrid = 0;
    }
    if (50 < AKTjetPt and AKTjetPt <= 100) {
        PTGrid = 1;
    }
    if (100 < AKTjetPt and AKTjetPt <= 150) {
        PTGrid = 2;
    }
    if (150 < AKTjetPt and AKTjetPt <= 200) {
        PTGrid = 3;
    }
    if (200 < AKTjetPt and AKTjetPt <= 250) {
        PTGrid = 4;
    }
    if (250 < AKTjetPt and AKTjetPt <= 300) {
        PTGrid = 5;
    }
    if (300 < AKTjetPt and AKTjetPt <= 350) {
        PTGrid = 6;
    }
    if (350 < AKTjetPt and AKTjetPt <= 400) {
        PTGrid = 7;
    }
    if (400 < AKTjetPt and AKTjetPt <= 450) {
        PTGrid = 8;
    }
    if (450 < AKTjetPt) {
        PTGrid = 9;
    }
    /*
    if (MuonTagging == true) {
        AKTjetPt = AKTjetPt/MuonwJER[ThetaGrid][PTGrid];
    } else {
        AKTjetPt = AKTjetPt/MuonwoJER[ThetaGrid][PTGrid];
    }
    */
    AKTjetPt = AKTjetPt / JER[ThetaGrid][PTGrid];
}

void Calibration(const char *inputFile,
    const char *outputFile,
        Double_t JER[10][10]) { 
    //Initiation
    gSystem -> Load("libDelphes.so");
    TFile *file_sig = new TFile(inputFile);
    TFile *output = new TFile(outputFile, "recreate");
    TTree *tree_sig = (TTree *) file_sig -> Get("Delphes");
    TTree *tree_output = new TTree("tree_output", "Delphes");

    cout << endl << "Initiating Calibration..." << endl;

    TLeaf *AKTjet_size = tree_sig -> GetLeaf("AKTjet_size");
    TLeaf *AKTjet_eta = tree_sig -> GetLeaf("AKTjet.Eta");
    TLeaf *AKTjet_phi = tree_sig -> GetLeaf("AKTjet.Phi");
    TLeaf *AKTjet_pt = tree_sig -> GetLeaf("AKTjet.PT");
    TLeaf *AKTjet_mass = tree_sig -> GetLeaf("AKTjet.Mass");

    TLeaf *GenJet_size = tree_sig -> GetLeaf("GenJet_size");
    TLeaf *GenJet_eta = tree_sig -> GetLeaf("GenJet.Eta");
    TLeaf *GenJet_phi = tree_sig -> GetLeaf("GenJet.Phi");
    TLeaf *GenJet_pt = tree_sig -> GetLeaf("GenJet.PT");
    TLeaf *GenJet_mass = tree_sig -> GetLeaf("GenJet.Mass");

    TLeaf *PID = tree_sig -> GetLeaf("Particle.PID");
    TLeaf *GenParticleEta = tree_sig -> GetLeaf("Particle.Eta");
    TLeaf *GenParticlePhi = tree_sig -> GetLeaf("Particle.Phi");
    TLeaf *GenParticlePt = tree_sig -> GetLeaf("Particle.PT");
    TLeaf *GenParticleSize = tree_sig -> GetLeaf("Particle_size");

    TLeaf *MuonEta = tree_sig -> GetLeaf("Muon.Eta");
    TLeaf *MuonPhi = tree_sig -> GetLeaf("Muon.Phi");
    TLeaf *MuonPt = tree_sig -> GetLeaf("Muon.PT");
    TLeaf *MuonSize = tree_sig -> GetLeaf("Muon_size");

    Int_t nEntries = tree_sig -> GetEntries();

    TH1D *AKTGenMass1Comp = new TH1D("AKTGenMass1Comp", "AKTGenMass1Comp", 50, -1, 2);
    TH2D *AKTjetPT_Theta = new TH2D("AKTjetPT_Theta", "AKTjetPT_Theta", 30, 0.1482, 3, 30, 0, 400);
    TH2D *GenJetPT_Theta = new TH2D("GenJetPT_Theta", "GenJetPT_Theta", 30, 0.1482, 3, 30, 0, 400);

    //Pairing variables
    Double_t AKTjet1eta1;
    Double_t AKTjet1theta1;
    Double_t AKTjet1phi1;
    Double_t AKTjet1pt1;
    Double_t AKTjet1mass1;
    Double_t AKTjet1eta2;
    Double_t AKTjet1theta2;
    Double_t AKTjet1phi2;
    Double_t AKTjet1pt2;
    Double_t AKTjet1mass2;
    Double_t AKTjet2eta1;
    Double_t AKTjet2theta1;
    Double_t AKTjet2phi1;
    Double_t AKTjet2pt1;
    Double_t AKTjet2mass1;
    Double_t AKTjet2eta2;
    Double_t AKTjet2theta2;
    Double_t AKTjet2phi2;
    Double_t AKTjet2pt2;
    Double_t AKTjet2mass2;
    Double_t AKTjetpairmass;

    Double_t Gen1eta;
    Double_t Gen1phi;
    Double_t Gen2eta;
    Double_t Gen2phi;
    Double_t Gen1pt;
    Double_t Gen2pt;
    Double_t Gen1mass;
    Double_t Gen2mass;

    Double_t Gen3eta;
    Double_t Gen3phi;
    Double_t Gen4eta;
    Double_t Gen4phi;
    Double_t Gen3pt;
    Double_t Gen4pt;
    Double_t Gen3mass;
    Double_t Gen4mass;

    Int_t pair1jet1entry;
    Int_t pair1jet2entry;
    Int_t pair2jet1entry;
    Int_t pair2jet2entry;

    Int_t AKTjet2entry1;
    Int_t AKTjet2entry2;

    Double_t diHiggsdis;

    Double_t AKTjetpair1Mass;
    Double_t AKTjetpair2Mass;

    Int_t AKTjetpaircnt;

    TLorentzVector AKTh1;
    TLorentzVector AKTh2;
    TLorentzVector AKTjet1;
    TLorentzVector AKTjet2;

    Double_t GenJetMass = 0;
    TLorentzVector Genh2;
    TLorentzVector GenJet1;
    TLorentzVector GenJet2;
    TLorentzVector GenJet3;
    TLorentzVector GenJet4;

    Double_t jet1DeltaR;
    Double_t jet2DeltaR;
    Double_t jet3DeltaR;
    Double_t jet4DeltaR;

    //Calibration variable
    Double_t AKTjetEta;
    Double_t AKTjetPt;
    Double_t AKTjetPhi;
    Double_t AKTjetMass;
    Double_t AKTjetTheta;

    Double_t GenJetEta;
    Double_t GenJetPt;
    Double_t GenJetPhi;
    Double_t GenJetTheta;

    Double_t GenMuonEta;
    Double_t GenMuonPhi;
    Double_t GenMuonPt;

    Double_t recoMuonEta;
    Double_t recoMuonPhi;
    Double_t recoMuonPt;

    TLorentzVector AKTjet;
    TLorentzVector GenJet;
    TLorentzVector recoMuon;
    TLorentzVector GenMuon;

    Double_t MatchedGenJetEta;
    Double_t MatchedGenJetPhi;

    Int_t ParticlePID;

    //Double_t JER[10][10];
    Double_t MuonwJER[10][10];
    Double_t MuonwoJER[10][10];

    Double_t matchJetCnt[10][10];
    Double_t matchJetMuonwCnt[10][10];
    Double_t matchJetMuonwoCnt[10][10];
    for (Int_t ThetaEntry = 0; ThetaEntry < 10; ThetaEntry++) {
        for (Int_t PTentry = 0; PTentry < 10; PTentry++) {
            matchJetCnt[ThetaEntry][PTentry] = 0;
            matchJetMuonwCnt[ThetaEntry][PTentry] = 0;
            matchJetMuonwoCnt[ThetaEntry][PTentry] = 0;
            JER[ThetaEntry][PTentry] = 0;
            MuonwJER[ThetaEntry][PTentry] = 0;
            MuonwoJER[ThetaEntry][PTentry] = 0;
        }
    }
    bool MuonTagging = false;

    Double_t JERtmp;
    Double_t DeltaPhi;
    Double_t DeltaR;
    Double_t MuonDeltaR;
    Double_t MuonDeltaPhi;
    Double_t GenMuonDeltaR;
    Double_t GenMuonDeltaPhi;
    Int_t PTGrid;
    Int_t ThetaGrid;

    //Calculate JER
    for (Long64_t entry = 0; entry < nEntries; entry++) {
        //Initiation
        tree_sig -> GetEntry(entry);
        AKTjet_size -> GetBranch() -> GetEntry(entry);
        GenJet_size -> GetBranch() -> GetEntry(entry);

        Int_t nAKTjet = AKTjet_size -> GetValue();
        Int_t nGenJet = GenJet_size -> GetValue();

        AKTjet_eta -> GetBranch() -> GetEntry(entry);
        AKTjet_phi -> GetBranch() -> GetEntry(entry);
        AKTjet_pt -> GetBranch() -> GetEntry(entry);
        AKTjet_mass -> GetBranch() -> GetEntry(entry);

        GenJet_eta -> GetBranch() -> GetEntry(entry);
        GenJet_pt -> GetBranch() -> GetEntry(entry);
        GenJet_mass -> GetBranch() -> GetEntry(entry);

        PID -> GetBranch() -> GetEntry(entry);
        GenParticleEta -> GetBranch() -> GetEntry(entry);
        GenParticlePhi -> GetBranch() -> GetEntry(entry);
        GenParticlePt -> GetBranch() -> GetEntry(entry);
        GenParticleSize -> GetBranch() -> GetEntry(entry);
        Int_t nParticle = GenParticleSize -> GetValue();

        MuonEta -> GetBranch() -> GetEntry(entry);
        MuonPhi -> GetBranch() -> GetEntry(entry);
        MuonPt -> GetBranch() -> GetEntry(entry);
        MuonSize -> GetBranch() -> GetEntry(entry);
        Int_t nMuon = MuonSize -> GetValue();
        if (nAKTjet == 6) {
            //Truth matching for all and calculate Jet energy response
            for (Int_t aktentry = 0; aktentry < nAKTjet; aktentry++) {
                AKTjetEta = AKTjet_eta -> GetValue(aktentry);
                AKTjetPhi = AKTjet_phi -> GetValue(aktentry);
                AKTjetPt = AKTjet_pt -> GetValue(aktentry);
                AKTjetMass = AKTjet_mass -> GetValue(aktentry);
                AKTjetTheta = 2 *atan(exp(-AKTjetEta));
                AKTjet.SetPtEtaPhiM(AKTjetPt, AKTjetEta, AKTjetPhi, AKTjetMass);
                DeltaR = 100;
                for (Int_t genentry = 0; genentry < nGenJet; genentry++) {
                    MuonTagging = false;
                    GenJetEta = GenJet_eta -> GetValue(genentry);
                    GenJetPhi = GenJet_phi -> GetValue(genentry);
                    GenJetPt = GenJet_pt -> GetValue(genentry);
                    //GenJetMass = GenJet_mass->GetValue(genentry)
                    GenJetTheta = 2 *atan(exp(-GenJetEta));

                    DeltaPhi = TMath::Abs(AKTjetPhi - GenJetPhi);
                    if (DeltaPhi > TMath::Pi()) {
                        DeltaPhi -= TMath::TwoPi();
                    }
                    Double_t DeltaRtmp = TMath::Sqrt(pow((AKTjetEta - GenJetEta), 2) + pow(DeltaPhi, 2));

                    //Double_t DeltaRtmp = GenJet.DeltaR(AKTjet);
                    if (DeltaRtmp < DeltaR) {
                        DeltaR = DeltaRtmp;
                        JERtmp = AKTjetPt / GenJetPt;
                        MatchedGenJetEta = GenJetEta;
                        MatchedGenJetPhi = GenJetPhi;
                        GenJet.SetPtEtaPhiM(GenJetPt, GenJetEta, GenJetPhi, 0);
                    }
                }

                if (DeltaR < 0.5) {
                    
                    if (0 <= AKTjetTheta and AKTjetTheta <= 0.3) {
                        ThetaGrid = 0;
                    }
                    if (0.3 < AKTjetTheta and AKTjetTheta <= 0.6) {
                        ThetaGrid = 1;
                    }
                    if (0.6 < AKTjetTheta and AKTjetTheta <= 0.9) {
                        ThetaGrid = 2;
                    }
                    if (0.9 < AKTjetTheta and AKTjetTheta <= 1.2) {
                        ThetaGrid = 3;
                    }
                    if (1.2 < AKTjetTheta and AKTjetTheta <= 1.5) {
                        ThetaGrid = 4;
                    }
                    if (1.5 < AKTjetTheta and AKTjetTheta <= 1.8) {
                        ThetaGrid = 5;
                    }
                    if (1.8 < AKTjetTheta and AKTjetTheta <= 2.1) {
                        ThetaGrid = 6;
                    }
                    if (2.1 < AKTjetTheta and AKTjetTheta <= 2.4) {
                        ThetaGrid = 7;
                    }
                    if (2.4 < AKTjetTheta and AKTjetTheta <= 2.7) {
                        ThetaGrid = 8;
                    }
                    if (2.7 < AKTjetTheta) {
                        ThetaGrid = 9;
                    }
                    if (0 <= AKTjetPt and AKTjetPt <= 50) {
                        PTGrid = 0;
                    }
                    if (50 < AKTjetPt and AKTjetPt <= 100) {
                        PTGrid = 1;
                    }
                    if (100 < AKTjetPt and AKTjetPt <= 150) {
                        PTGrid = 2;
                    }
                    if (150 < AKTjetPt and AKTjetPt <= 200) {
                        PTGrid = 3;
                    }
                    if (200 < AKTjetPt and AKTjetPt <= 250) {
                        PTGrid = 4;
                    }
                    if (250 < AKTjetPt and AKTjetPt <= 300) {
                        PTGrid = 5;
                    }
                    if (300 < AKTjetPt and AKTjetPt <= 350) {
                        PTGrid = 6;
                    }
                    if (350 < AKTjetPt and AKTjetPt <= 400) {
                        PTGrid = 7;
                    }
                    if (400 < AKTjetPt and AKTjetPt <= 450) {
                        PTGrid = 8;
                    }
                    if (450 < AKTjetPt) {
                        PTGrid = 9;
                    }
                    matchJetCnt[ThetaGrid][PTGrid] += 1;
                    JER[ThetaGrid][PTGrid] = (JER[ThetaGrid][PTGrid] *(matchJetCnt[ThetaGrid][PTGrid] - 1) + JERtmp) / matchJetCnt[ThetaGrid][PTGrid];
		    
                } else {
                    //cout << DeltaR << endl;
                }
            }
        }
    }
    cout << endl << "Calibration's result:" << endl;
    //Calibration result
    TCanvas *myCaliCanvas = new TCanvas("myCaliCanvas", "My Canvas", 800, 600);
    //myCaliCanvas -> SaveAs("jetPTresponse.png");
    
    for (Int_t PTentry = 9; PTentry >= 0; PTentry -= 1) {
        for (Int_t ThetaEntry = 0; ThetaEntry < 10; ThetaEntry++) {
            if (matchJetCnt[ThetaEntry][PTentry] == 0) {
                JER[ThetaEntry][PTentry] = 1.0000;
            }
        }
    }
    cout << "Jet Energy Calibration gives following Jet Energy Response:" << endl;
    for (Int_t PTentry = 9; PTentry >= 0; PTentry -= 1) {
        for (Int_t ThetaEntry = 0; ThetaEntry < 10; ThetaEntry++) {
            cout << floorf(JER[ThetaEntry][PTentry] *100) / 100 << ", ";
        }
        cout << endl;
    }
    
    cout << endl << "Number of Jets included in calculation" << endl;
    for (Int_t PTentry = 9; PTentry >= 0; PTentry -= 1) {
        for (Int_t ThetaEntry = 0; ThetaEntry < 10; ThetaEntry++) {
            cout << matchJetCnt[ThetaEntry][PTentry] << ", ";
        }
        cout << endl;
    }
    

    tree_output -> Write();
    output -> Close();
    file_sig -> Close();

    cout << endl << "Calibration all done! Quit..." << endl;
}

void Pairing_w_JES(const char *inputFile,
    Double_t JER[10][10],
    	TH1D* AKTjetMass1,
	        TH1D* AKTjetMass2,
		    TH2F* AKTjetPairs,
	    	        TTree* tree_BDT,
		            TFile* file_BDT) {
    //clock_t tStart = clock();
    //Initiation
    gSystem -> Load("libDelphes.so");
    gROOT -> ProcessLine("gErrorIgnoreLevel = kWarning;");
    gErrorIgnoreLevel = kWarning;
    //TChain chain("Delphes");
    //chain.Add(inputFile);
    //ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
    TFile *file_sig = new TFile(inputFile);
    //TFile *output = new TFile(outputFile, "RECREATE");
    TTree *tree_sig = (TTree*) file_sig -> Get("Delphes");
    //TTree *tree_output = new TTree("tree_output", "Delphes");

    cout << endl << "///////////////////////////////////////////////////" << endl << "Initiating pairing algorithm..." << endl;

    TLeaf *Lep_mu_size = tree_sig -> GetLeaf("Muon_size");
    TLeaf *Lep_mu_eta = tree_sig -> GetLeaf("Muon.Eta");
    TLeaf *Lep_mu_phi = tree_sig -> GetLeaf("Muon.Phi");
    TLeaf *Lep_mu_pt = tree_sig -> GetLeaf("Muon.PT");
    TLeaf *Lep_mu_Charge = tree_sig -> GetLeaf("Muon.Charge");

    TLeaf *Lep_e_size = tree_sig -> GetLeaf("Electron_size");
    TLeaf *Lep_e_eta = tree_sig -> GetLeaf("Electron.Eta");
    TLeaf *Lep_e_phi = tree_sig -> GetLeaf("Electron.Phi");
    TLeaf *Lep_e_pt = tree_sig -> GetLeaf("Electron.PT");
    TLeaf *Lep_e_Charge = tree_sig -> GetLeaf("Electron.Charge");

    TLeaf *AKTjet_size = tree_sig -> GetLeaf("AKTjet_size");
    TLeaf *AKTjet_eta = tree_sig -> GetLeaf("AKTjet.Eta");
    TLeaf *AKTjet_phi = tree_sig -> GetLeaf("AKTjet.Phi");
    TLeaf *AKTjet_pt = tree_sig -> GetLeaf("AKTjet.PT");
    TLeaf *AKTjet_mass = tree_sig -> GetLeaf("AKTjet.Mass");
    TLeaf *AKTjet_BTag = tree_sig -> GetLeaf("AKTjet.BTag");
    TLeaf *AKTjet_TauTag = tree_sig -> GetLeaf("AKTjet.TauTag");
    TLeaf *AKTjet_Charge = tree_sig -> GetLeaf("AKTjet.Charge");

    TLeaf *GenJet_size = tree_sig -> GetLeaf("GenJet_size");
    TLeaf *GenJet_eta = tree_sig -> GetLeaf("GenJet.Eta");
    TLeaf *GenJet_phi = tree_sig -> GetLeaf("GenJet.Phi");
    TLeaf *GenJet_pt = tree_sig -> GetLeaf("GenJet.PT");
    TLeaf *GenJet_mass = tree_sig -> GetLeaf("GenJet.Mass");

    TLeaf *PID = tree_sig -> GetLeaf("Particle.PID");
    TLeaf *GenParticleEta = tree_sig -> GetLeaf("Particle.Eta");
    TLeaf *GenParticlePhi = tree_sig -> GetLeaf("Particle.Phi");
    TLeaf *GenParticlePt = tree_sig -> GetLeaf("Particle.PT");
    TLeaf *GenParticleSize = tree_sig -> GetLeaf("Particle_size");

    TLeaf *MissingETMET = tree_sig -> GetLeaf("MissingET.MET");
    TLeaf *MissingETEta = tree_sig -> GetLeaf("MissingET.Eta");
    TLeaf *MissingETPhi = tree_sig -> GetLeaf("MissingET.Phi");

    Int_t nEntries = tree_sig -> GetEntries();

    //TH1D *AKTjetMass1 = new TH1D("AKTjetMass1", "Anti_KTjet leading jets pair invariant mass", 150, 0, 400);
    //TH1D *AKTjetMass2 = new TH1D("AKTjetMass2", "Anti_KTjet sub-leading jets pair invariant mass", 150, 0, 400);
    TH1D *GenAKTMass2 = new TH1D("GenAKTMass2", "GenAKTMass2", 50, 0, 600);
    TH1D *AKTGenMass1Comp = new TH1D("AKTGenMass1Comp", "AKTGenMass1Comp", 50, -1, 2);
    TH1D *AKTGenMass2Comp = new TH1D("AKTGenMass2Comp", "AKTGenMass2Comp", 50, -1, 2);
    TH1D *GenUncutMass2 = new TH1D("GenUncutMass2", "GenUncutMass2", 50, 0, 600);

    TH1D *AKTGenPt1Comp = new TH1D("AKTGenPt1Comp", "AKTGenPt1Comp", 200, -1, 5);
    TH1D *AKTGenPt2Comp = new TH1D("AKTGenPt2Comp", "AKTGenPt2Comp", 200, -1, 5);
    TH1D *AKTGenPt3Comp = new TH1D("AKTGenPt3Comp", "AKTGenPt3Comp", 200, -1, 5);
    TH1D *AKTGenPt4Comp = new TH1D("AKTGenPt4Comp", "AKTGenPt4Comp", 200, -1, 5);

    TH2D *AKTjetPT_Theta = new TH2D("AKTjetPT_Theta", "AKTjetPT_Theta", 30, 0.1482, 3, 30, 0, 400);
    TH2D *GenJetPT_Theta = new TH2D("GenJetPT_Theta", "GenJetPT_Theta", 30, 0.1482, 3, 30, 0, 400);

    Double_t m_e;
    m_e = 1.777;
    Double_t m_mu;
    m_mu = 1.777;
    Double_t LEPeta;
    Double_t LEPphi;
    Double_t LEPpt;
    Double_t LEPmass;
    Int_t LEPcharge;

    Double_t AKTjet1eta2;
    Double_t AKTjet1theta2;
    Double_t AKTjet1phi2;
    Double_t AKTjet1pt2;
    Double_t AKTjet1mass2;
    Int_t AKTjet1charge2;

    Double_t AKTjet2eta1;
    Double_t AKTjet2theta1;
    Double_t AKTjet2phi1;
    Double_t AKTjet2pt1;
    Double_t AKTjet2mass1;
    
    Double_t AKTjet2eta2;
    Double_t AKTjet2theta2;
    Double_t AKTjet2phi2;
    Double_t AKTjet2pt2;
    Double_t AKTjet2mass2;
    Double_t AKTjetpairmass;

    Int_t AKTpairCharge;

    Int_t AKTjet1BTag2;
    Int_t AKTjet2BTag1;
    Int_t AKTjet2BTag2;

    Int_t AKTjet1TauTag2;
    Int_t AKTjet2TauTag1;
    Int_t AKTjet2TauTag2;

    Double_t Tau1eta;
    Double_t Tau1phi;
    Double_t Tau1pt;
    Double_t Tau2eta;
    Double_t Tau2phi;
    Double_t Tau2pt;

    Int_t Epair1jet1entry;
    Int_t Epair1jet2entry;
    Int_t Mupair1jet1entry;
    Int_t Mupair1jet2entry;
    Int_t pair1jet1entry;
    Int_t pair1jet2entry;
    Int_t pair2jet1entry;
    Int_t pair2jet2entry;
/*
    Int_t AKTjet2entry1;
    Int_t AKTjet2entry2;
*/
    Double_t diHiggsdis;

    Double_t Ejetpair1Mass;
    Double_t Mujetpair1Mass;

    Double_t AKTjetpair1Mass;
    Double_t AKTjetpair2Mass;

    Int_t EAKTjetpaircnt;
    Int_t MuAKTjetpaircnt;

    TLorentzVector AKTh1;
    TLorentzVector AKTh2;
    TLorentzVector AKT1jet1;
    TLorentzVector AKT1jet2;
    TLorentzVector AKT2jet1;
    TLorentzVector AKT2jet2;

    //Calibration variable
    Double_t AKTjetEta;
    Double_t AKTjetPt;
    Double_t AKTjetPhi;
    Double_t AKTjetMass;
    Double_t AKTjetTheta;

    Double_t GenJetEta;
    Double_t GenJetPt;
    Double_t GenJetPhi;
    Double_t GenJetTheta;

    Double_t GenMuonEta;
    Double_t GenMuonPhi;
    Double_t GenMuonPt;

    Double_t recoMuonEta;
    Double_t recoMuonPhi;
    Double_t recoMuonPt;

    TLorentzVector AKTjet;
    TLorentzVector GenJet;
    TLorentzVector recoMuon;
    TLorentzVector GenMuon;

    Double_t MatchedGenJetEta;
    Double_t MatchedGenJetPhi;

    Int_t ParticlePID;

    //BDT variable
    Int_t BDTNjets;

    Double_t BDTlep1pt1;
    Double_t BDTjet1pt2;
    Double_t BDTjet2pt1;
    Double_t BDTjet2pt2;

    Double_t BDTlep1eta1;
    Double_t BDTjet1eta2;
    Double_t BDTjet2eta1;
    Double_t BDTjet2eta2;

    Double_t BDTlep1phi1;
    Double_t BDTjet1phi2;
    Double_t BDTjet2phi1;
    Double_t BDTjet2phi2;

    Int_t BDTlep1charge1;
    Int_t BDTjet1charge2;

    Int_t BDTpairCharge;

    Int_t BDTjet1BTag2;
    Int_t BDTjet2BTag1;
    Int_t BDTjet2BTag2;

    Int_t BDTjet1TauTag2;
    Int_t BDTjet2TauTag1;
    Int_t BDTjet2TauTag2;

    Double_t BDThiggs1pt;
    Double_t BDThiggs1eta;
    Double_t BDThiggs1phi;
    Double_t BDThiggs1invm;

    Double_t BDThiggs2pt;
    Double_t BDThiggs2eta;
    Double_t BDThiggs2phi;
    Double_t BDThiggs2invm;

    Double_t BDThiggsdeltaEta;
    Double_t BDThiggsdeltaPhi;
    
    Double_t BDTMissingETMET;
    Double_t BDTMissingETEta;
    Double_t BDTMissingETPhi;

    Double_t BDTdihiggspt;
    Double_t BDTdihiggseta;
    Double_t BDTdihiggsphi;
    Double_t BDTdihiggsinvm;

    Double_t BDThiggsAngle;
    
    Double_t DiHiggsMETDeltaTheta;

    Int_t NBTag;
    Int_t NTauTag;

    tree_BDT -> Branch("BDTNjets", &BDTNjets, "BDTNjets/I");
    tree_BDT -> Branch("NBTag", &NBTag, "NBTag/I");
    tree_BDT -> Branch("NTauTag", &NTauTag, "NTauTag/I");

    tree_BDT -> Branch("BDTjet1BTag2", &BDTjet1BTag2, "BDTjet1BTag2/I");
    tree_BDT -> Branch("BDTjet2BTag1", &BDTjet2BTag1, "BDTjet2BTag1/I");
    tree_BDT -> Branch("BDTjet2BTag2", &BDTjet2BTag2, "BDTjet2BTag2/I");

    tree_BDT -> Branch("BDTjet1TauTag2", &BDTjet1TauTag2, "BDTjet1TauTag2/I");
    tree_BDT -> Branch("BDTjet2TauTag1", &BDTjet2TauTag1, "BDTjet2TauTag1/I");
    tree_BDT -> Branch("BDTjet2TauTag2", &BDTjet2TauTag2, "BDTjet2TauTag2/I");

    tree_BDT -> Branch("BDTlep1pt1", &BDTlep1pt1, "BDTlep1pt1/D");
    tree_BDT -> Branch("BDTjet1pt2", &BDTjet1pt2, "BDTjet1pt2/D");
    tree_BDT -> Branch("BDTjet2pt1", &BDTjet2pt1, "BDTjet2pt1/D");
    tree_BDT -> Branch("BDTjet2pt2", &BDTjet2pt2, "BDTjet2pt2/D");

    tree_BDT -> Branch("BDTlep1eta1", &BDTlep1eta1, "BDTlep1eta1/D");
    tree_BDT -> Branch("BDTjet1eta2", &BDTjet1eta2, "BDTjet1eta2/D");
    tree_BDT -> Branch("BDTjet2eta1", &BDTjet2eta1, "BDTjet2eta1/D");
    tree_BDT -> Branch("BDTjet2eta2", &BDTjet2eta2, "BDTjet2eta2/D");

    tree_BDT -> Branch("BDTlep1phi1", &BDTlep1phi1, "BDTlep1phi1/D");
    tree_BDT -> Branch("BDTjet1phi2", &BDTjet1phi2, "BDTjet1phi2/D");
    tree_BDT -> Branch("BDTjet2phi1", &BDTjet2phi1, "BDTjet2phi1/D");
    tree_BDT -> Branch("BDTjet2phi2", &BDTjet2phi2, "BDTjet2phi2/D");

    tree_BDT -> Branch("BDTlep1charge1", &BDTlep1charge1, "BDTlep1charge1/I");
    tree_BDT -> Branch("BDTjet1charge2", &BDTjet1charge2, "BDTjet1charge2/I");

    tree_BDT -> Branch("BDTpairCharge", &BDTpairCharge, "BDTpairCharge/I");

    tree_BDT -> Branch("BDThiggs1pt", &BDThiggs1pt, "BDThiggs1pt/D");
    tree_BDT -> Branch("BDThiggs1eta", &BDThiggs1eta, "BDThiggs1eta/D");
    tree_BDT -> Branch("BDThiggs1phi", &BDThiggs1phi, "BDThiggs1phi/D");
    tree_BDT -> Branch("BDThiggs1invm", &BDThiggs1invm, "BDThiggs1invm/D");

    tree_BDT -> Branch("BDThiggs2pt", &BDThiggs2pt, "BDThiggs2pt/D");
    tree_BDT -> Branch("BDThiggs2eta", &BDThiggs2eta, "BDThiggs2eta/D");
    tree_BDT -> Branch("BDThiggs2phi", &BDThiggs2phi, "BDThiggs2phi/D");
    tree_BDT -> Branch("BDThiggs2invm", &BDThiggs2invm, "BDThiggs2invm/D");

    tree_BDT -> Branch("BDThiggsdeltaEta", &BDThiggsdeltaEta, "BDThiggsdeltaEta/D");
    tree_BDT -> Branch("BDThiggsdeltaPhi", &BDThiggsdeltaPhi, "BDThiggsdeltaPhi/D");
    
    tree_BDT -> Branch("BDThiggsAngle", &BDThiggsAngle, "BDThiggsAngle/D");
    
    tree_BDT -> Branch("DiHiggsMETDeltaTheta", &DiHiggsMETDeltaTheta, "DiHiggsMETDeltaTheta/D");

    tree_BDT -> Branch("BDTMissingETMET", &BDTMissingETMET, "BDTMissingETMET/D");
    tree_BDT -> Branch("BDTMissingETEta", &BDTMissingETEta, "BDTMissingETEta/D");
    tree_BDT -> Branch("BDTMissingETPhi", &BDTMissingETPhi, "BDTMissingETPhi/D");

    tree_BDT -> Branch("BDTdihiggspt", &BDTdihiggspt, "BDTdihiggspt/D");
    tree_BDT -> Branch("BDTdihiggseta", &BDTdihiggseta, "BDTdihiggseta/D");
    tree_BDT -> Branch("BDTdihiggsphi", &BDTdihiggsphi, "BDTdihiggsphi/D");
    tree_BDT -> Branch("BDTdihiggsinvm", &BDTdihiggsinvm, "BDTdihiggsinvm/D");

    //cout << endl << "Calling Calibration algorithm..." << endl;
    //Calibration(inputFileForJES, outputFileForJES, JER); // MuonwJER, MuonwoJER);

    //run over event entries
    for (Long64_t entry = 0; entry < nEntries; entry++) {
        //Initiation

	//cout << "Entry = " << entry << endl;
        tree_sig -> GetEntry(entry);
        //tree_output -> GetEntry(entry);
        AKTjet_size -> GetBranch() -> GetEntry(entry);
        Lep_e_size -> GetBranch() -> GetEntry(entry);
        Lep_mu_size -> GetBranch() -> GetEntry(entry);
        GenParticleSize -> GetBranch() -> GetEntry(entry);

        Int_t nAKTjet = AKTjet_size -> GetValue();
        Int_t nElectron = Lep_e_size -> GetValue();
        Int_t nMuon = Lep_mu_size -> GetValue();
        Int_t nParticle = GenParticleSize -> GetValue();

	if (nMuon == 0 and nElectron == 0) {
	    continue;
	}
	//cout << "Event # = " << entry << "; AKTjet = " << nAKTjet << "; nElectron = " << nElectron << "; nMuon = " << nMuon << endl;

        AKTjet_eta -> GetBranch() -> GetEntry(entry);
        AKTjet_phi -> GetBranch() -> GetEntry(entry);
        AKTjet_pt -> GetBranch() -> GetEntry(entry);
        AKTjet_mass -> GetBranch() -> GetEntry(entry);
        AKTjet_BTag -> GetBranch() -> GetEntry(entry);
	AKTjet_TauTag -> GetBranch() -> GetEntry(entry);
	AKTjet_Charge -> GetBranch() -> GetEntry(entry);

	PID -> GetBranch() -> GetEntry(entry);
	GenParticleEta -> GetBranch() -> GetEntry(entry);
	GenParticlePhi -> GetBranch() -> GetEntry(entry);
	GenParticlePt -> GetBranch() -> GetEntry(entry);

	Lep_e_eta -> GetBranch() -> GetEntry(entry);
	Lep_e_phi -> GetBranch() -> GetEntry(entry);
	Lep_e_pt -> GetBranch() -> GetEntry(entry);
	Lep_e_Charge -> GetBranch() -> GetEntry(entry);

	Lep_mu_eta -> GetBranch() -> GetEntry(entry);
	Lep_mu_phi -> GetBranch() -> GetEntry(entry);
	Lep_mu_pt -> GetBranch() -> GetEntry(entry);
	Lep_mu_Charge -> GetBranch() -> GetEntry(entry);

	MissingETMET -> GetBranch() -> GetEntry(entry);
	MissingETEta -> GetBranch() -> GetEntry(entry);
	MissingETPhi -> GetBranch() -> GetEntry(entry);

	BDTMissingETMET = MissingETMET -> GetValue();
        BDTMissingETEta = MissingETEta -> GetValue();
        BDTMissingETPhi = MissingETPhi -> GetValue();

        AKTjetpair1Mass = 10000;
        AKTjetpair2Mass = 10000;

	Ejetpair1Mass = 0;
	Mujetpair1Mass = 0;
        Epair1jet1entry = 0;
        Epair1jet2entry = 0;
        Mupair1jet1entry = 0;
        Mupair1jet2entry = 0;
        pair1jet1entry = 0;
        pair1jet2entry = 0;
        pair2jet1entry = 0;
        pair2jet2entry = 0;

        diHiggsdis = 1000000000;

        BDTjet1BTag2 = 0;
        BDTjet2BTag1 = 0;
        BDTjet2BTag2 = 0;
        
        BDTjet1TauTag2 = 0;
        BDTjet2TauTag1 = 0;
        BDTjet2TauTag2 = 0;

        //Pairing up Anti-kt jets
        //Double_t AKTjetpair[nAKTjet * (nAKTjet - 1) / 2][3];
	Double_t EJetPair[nAKTjet * nElectron][3];
        Double_t MuJetPair[nAKTjet * nMuon][3];

        EAKTjetpaircnt = 0;
        MuAKTjetpaircnt = 0;

        if (nAKTjet >= 3 and (nElectron + nMuon) > 0) {
            /*
	    for (Int_t aktentry = 0; aktentry < nAKTjet; aktentry++) {
                AKTjetEta = AKTjet_eta->GetValue(aktentry);
                AKTjetPt = AKTjet_pt->GetValue(aktentry);
                AKTjetTheta = 2 * atan(exp(-AKTjetEta));
                AKTjetPT_Theta->Fill(AKTjetTheta, AKTjetPt);
                JetEnergyFix(AKTjetEta, AKTjetPt, JER);
		
            }
    */

  	    TVector3 bjet1;
	    TVector3 bjet2;
	    AKTjetpair2Mass = 0;
            for (Int_t akt1entry2 = 0; akt1entry2 < nAKTjet; akt1entry2++) {
		AKTjet2BTag1 = AKTjet_BTag -> GetValue(akt1entry2);
                if (/*akt1entry2 != pair1jet2entry and*/ AKTjet2BTag1 >= 4) {
                    for (Int_t akt2entry2 = 0; akt2entry2 < nAKTjet; akt2entry2++) {
			AKTjet2BTag2 = AKTjet_BTag -> GetValue(akt2entry2);
/*
			AKTjet2eta1 = AKTjet_eta -> GetValue(pair2jet1entry);
			AKTjet2phi1 = AKTjet_phi -> GetValue(pair2jet1entry);
			AKTjet2eta2 = AKTjet_eta -> GetValue(pair2jet2entry);
			AKTjet2phi2 = AKTjet_phi -> GetValue(pair2jet2entry);
			bjet1.SetPtEtaPhi(1.0, AKTjet2eta1, AKTjet2phi1);

			bjet2.SetPtEtaPhi(1.0, AKTjet2eta2, AKTjet2phi2);
                        */
			if (akt2entry2 != akt1entry2 and AKTjet2BTag2 >= 4) {

                            AKTjet2eta1 = AKTjet_eta -> GetValue(akt1entry2);
                            AKTjet2phi1 = AKTjet_phi -> GetValue(akt1entry2);
                            AKTjet2pt1 = AKTjet_pt -> GetValue(akt1entry2);
                            JetEnergyFix(AKTjet2eta1, AKTjet2pt1, JER);
                            AKTjet2mass1 = AKTjet_mass -> GetValue(akt1entry2);
		
                            AKTjet2eta2 = AKTjet_eta -> GetValue(akt2entry2);
                            AKTjet2phi2 = AKTjet_phi -> GetValue(akt2entry2);
                            AKTjet2pt2 = AKTjet_pt -> GetValue(akt2entry2);
                            JetEnergyFix(AKTjet2eta2, AKTjet2pt2, JER);
                            AKTjet2mass2 = AKTjet_mass -> GetValue(akt2entry2);
	
                            AKT2jet1.SetPtEtaPhiM(AKTjet2pt1, AKTjet2eta1, AKTjet2phi1, AKTjet2mass1);
                            AKT2jet2.SetPtEtaPhiM(AKTjet2pt2, AKTjet2eta2, AKTjet2phi2, AKTjet2mass2);
                            AKTh2 = AKT2jet1 + AKT2jet2;
                            AKTjetpairmass = AKTh2.Mag();

                            if (TMath::Abs(125 - AKTjetpairmass) < TMath::Abs(125 - AKTjetpair2Mass)) {
                                AKTjetpair2Mass = AKTjetpairmass;
                                pair2jet1entry = akt1entry2;
                                pair2jet2entry = akt2entry2;

                            }
                        }
                    }
                }
            }

	    if (pair2jet1entry == 0 and pair2jet2entry == 0) {
	        continue;
	    }
	    AKTjet2eta1 = AKTjet_eta -> GetValue(pair2jet1entry);
	    AKTjet2phi1 = AKTjet_phi -> GetValue(pair2jet1entry);
	    AKTjet2eta2 = AKTjet_eta -> GetValue(pair2jet2entry);
	    AKTjet2phi2 = AKTjet_phi -> GetValue(pair2jet2entry);

	    bjet1.SetPtEtaPhi(1.0, AKTjet2eta1, AKTjet2phi1);
	    bjet2.SetPtEtaPhi(1.0, AKTjet2eta2, AKTjet2phi2);

            //Pairing up leading Tau_Tagged anti-KT jet pair
            for (Int_t akt1entry = 0; akt1entry < nElectron; akt1entry++) {
                for (Int_t akt2entry = 0; akt2entry < nAKTjet; akt2entry++) {
		    
                    AKTjet1TauTag2 = AKTjet_TauTag -> GetValue(akt2entry);
                    AKTjet1BTag2 = AKTjet_BTag -> GetValue(akt2entry);
		    AKTjet1eta2 = AKTjet_eta -> GetValue(akt2entry);
                    AKTjet1phi2 = AKTjet_phi -> GetValue(akt2entry);

		    TVector3 jet2;
		    TVector3 lep2;

		    LEPcharge = Lep_e_Charge -> GetValue(akt1entry);
		    AKTjet1charge2 = AKTjet_Charge -> GetValue(akt2entry);
		
		    AKTpairCharge = LEPcharge * AKTjet1charge2;

		    LEPeta = Lep_e_eta -> GetValue(akt1entry);
		    LEPphi = Lep_e_phi -> GetValue(akt1entry);
		    
		    jet2.SetPtEtaPhi(1.0, AKTjet1eta2, AKTjet1phi2);
		    lep2.SetPtEtaPhi(1.0, LEPeta, LEPphi);
                    if (AKTjet1TauTag2 == 1 and AKTpairCharge == -1 and /* lep2.DeltaR(bjet1) > 1.0 and lep2.DeltaR(bjet2) > 1.0 and*/ AKTjet1BTag2 == 0/* and akt2entry != pair2jet1entry and akt2entry != pair2jet2entry*/) {

                        LEPeta = Lep_e_eta -> GetValue(akt1entry);
                        LEPphi = Lep_e_phi -> GetValue(akt1entry);
                        LEPpt = Lep_e_pt -> GetValue(akt1entry);
                        LEPcharge = Lep_e_Charge -> GetValue(akt1entry);
		
                        AKTjet1eta2 = AKTjet_eta -> GetValue(akt2entry);
                        AKTjet1phi2 = AKTjet_phi -> GetValue(akt2entry);
                        AKTjet1pt2 = AKTjet_pt -> GetValue(akt2entry);
                        JetEnergyFix(AKTjet1eta2, AKTjet1pt2, JER);
                        AKTjet1mass2 = AKTjet_mass -> GetValue(akt2entry);
                        AKTjet1charge2 = AKTjet_Charge -> GetValue(akt2entry);
		
                        AKTjetpairmass = 0;

                        AKT1jet1.SetPtEtaPhiM(LEPpt, LEPeta, LEPphi, m_e);
                        AKT1jet2.SetPtEtaPhiM(AKTjet1pt2, AKTjet1eta2, AKTjet1phi2, AKTjet1mass2);
                        AKTh1 = AKT1jet1 + AKT1jet2;
                        AKTjetpairmass = AKTh1.Mag();

                        EJetPair[EAKTjetpaircnt][0] = AKTjetpairmass;
                        EJetPair[EAKTjetpaircnt][1] = akt1entry;
                        EJetPair[EAKTjetpaircnt][2] = akt2entry;
/*
			if (EJetPair[EAKTjetpaircnt][2] > 10) {
			    cout << "E akt eird entry number! "<< akt2entry <<endl;
			}
*/
                    }
                    EAKTjetpaircnt++;
                }
            }
            for (Int_t akt1entry = 0; akt1entry < nMuon; akt1entry++) {
                for (Int_t akt2entry = 0; akt2entry < nAKTjet; akt2entry++) {

                    AKTjet1TauTag2 = AKTjet_TauTag -> GetValue(akt2entry);
                    AKTjet1BTag2 = AKTjet_BTag -> GetValue(akt2entry);
                    AKTjet1eta2 = AKTjet_eta -> GetValue(akt2entry);
                    AKTjet1phi2 = AKTjet_phi -> GetValue(akt2entry);

		    TVector3 jet1;
		    TVector3 lep1;

                    LEPcharge = Lep_mu_Charge -> GetValue(akt1entry);
                    AKTjet1charge2 = AKTjet_Charge -> GetValue(akt2entry);
                    AKTpairCharge = LEPcharge * AKTjet1charge2;
		
		    LEPeta = Lep_mu_eta -> GetValue(akt1entry);
		    LEPphi = Lep_mu_phi -> GetValue(akt1entry);
		    
		    jet1.SetPtEtaPhi(1.0, AKTjet1eta2, AKTjet1phi2);
		    lep1.SetPtEtaPhi(1.0, LEPeta, LEPphi);
                    if (AKTjet1TauTag2 == 1 and AKTpairCharge == -1 /*and lep1.DeltaR(bjet1) > 1.0 and lep1.DeltaR(bjet2) > 1.0*/ and AKTjet1BTag2 == 0/* and akt2entry != pair2jet1entry and akt2entry != pair2jet2entry*/) {

                        LEPeta = Lep_mu_eta -> GetValue(akt1entry);
                        LEPphi = Lep_mu_phi -> GetValue(akt1entry);
                        LEPpt = Lep_mu_pt -> GetValue(akt1entry);
                        LEPcharge = Lep_mu_Charge -> GetValue(akt1entry);
		
                        AKTjet1eta2 = AKTjet_eta -> GetValue(akt2entry);
                        AKTjet1phi2 = AKTjet_phi -> GetValue(akt2entry);
                        AKTjet1pt2 = AKTjet_pt -> GetValue(akt2entry);
                        JetEnergyFix(AKTjet1eta2, AKTjet1pt2, JER);
                        AKTjet1mass2 = AKTjet_mass -> GetValue(akt2entry);
                        AKTjet1charge2 = AKTjet_Charge -> GetValue(akt2entry);
	
                        AKTjetpairmass = 0;

                        AKT1jet1.SetPtEtaPhiM(LEPpt, LEPeta, LEPphi, m_mu);
                        AKT1jet2.SetPtEtaPhiM(AKTjet1pt2, AKTjet1eta2, AKTjet1phi2, AKTjet1mass2);
                        AKTh1 = AKT1jet1 + AKT1jet2;
                        AKTjetpairmass = AKTh1.Mag();

                        MuJetPair[MuAKTjetpaircnt][0] = AKTjetpairmass;
                        MuJetPair[MuAKTjetpaircnt][1] = akt1entry;
                        MuJetPair[MuAKTjetpaircnt][2] = akt2entry;

                    }
                    MuAKTjetpaircnt++;
                }
            }

            Ejetpair1Mass = 0;
            for (Int_t pairindex = 0; pairindex < (nAKTjet * nElectron); pairindex++) {
                if (TMath::Abs(EJetPair[pairindex][0] - 125) < TMath::Abs(Ejetpair1Mass - 125)) {
                    Ejetpair1Mass = EJetPair[pairindex][0];
                    Epair1jet1entry = EJetPair[pairindex][1];
                    Epair1jet2entry = EJetPair[pairindex][2];

                }
            }
            Mujetpair1Mass = 0;
            for (Int_t pairindex = 0; pairindex < (nAKTjet * nMuon); pairindex++) {
                if (TMath::Abs(MuJetPair[pairindex][0] - 125) < TMath::Abs(Mujetpair1Mass - 125)) {
                    Mujetpair1Mass = MuJetPair[pairindex][0];
                    Mupair1jet1entry = MuJetPair[pairindex][1];
                    Mupair1jet2entry = MuJetPair[pairindex][2];

                }
            }

	    Int_t LepFlag = 2;

            if (TMath::Abs(Ejetpair1Mass - 125) < TMath::Abs(Mujetpair1Mass - 125)) {
                LepFlag = 0;
                AKTjetpair1Mass = Ejetpair1Mass;
                pair1jet1entry = Epair1jet1entry;
                pair1jet2entry = Epair1jet2entry;
            } else if (TMath::Abs(Ejetpair1Mass - 125) >= TMath::Abs(Mujetpair1Mass - 125) and Mujetpair1Mass > 0) {
                LepFlag = 1;
                AKTjetpair1Mass = Mujetpair1Mass;
                pair1jet1entry = Mupair1jet1entry;
                pair1jet2entry = Mupair1jet2entry;
	    } else {
		continue;
	    }

	    /*
            AKTjetpair2Mass = 0;
            for (Int_t akt1entry2 = 0; akt1entry2 < nAKTjet; akt1entry2++) {
		AKTjet2BTag1 = AKTjet_BTag -> GetValue(akt1entry2);
                if (akt1entry2 != pair1jet2entry and AKTjet2BTag1 >= 4) {
                    for (Int_t akt2entry2 = 0; akt2entry2 < nAKTjet; akt2entry2++) {
			AKTjet2BTag2 = AKTjet_BTag -> GetValue(akt2entry2);
                        if (akt2entry2 != pair1jet2entry and akt2entry2 != akt1entry2 and AKTjet2BTag2 >= 4) {

                            AKTjet2eta1 = AKTjet_eta -> GetValue(akt1entry2);
                            AKTjet2phi1 = AKTjet_phi -> GetValue(akt1entry2);
                            AKTjet2pt1 = AKTjet_pt -> GetValue(akt1entry2);
                            JetEnergyFix(AKTjet2eta1, AKTjet2pt1, JER);
                            AKTjet2mass1 = AKTjet_mass -> GetValue(akt1entry2);
		
                            AKTjet2eta2 = AKTjet_eta -> GetValue(akt2entry2);
                            AKTjet2phi2 = AKTjet_phi -> GetValue(akt2entry2);
                            AKTjet2pt2 = AKTjet_pt -> GetValue(akt2entry2);
                            JetEnergyFix(AKTjet2eta2, AKTjet2pt2, JER);
                            AKTjet2mass2 = AKTjet_mass -> GetValue(akt2entry2);
	
                            AKT2jet1.SetPtEtaPhiM(AKTjet2pt1, AKTjet2eta1, AKTjet2phi1, AKTjet2mass1);
                            AKT2jet2.SetPtEtaPhiM(AKTjet2pt2, AKTjet2eta2, AKTjet2phi2, AKTjet2mass2);
                            AKTh2 = AKT2jet1 + AKT2jet2;
                            AKTjetpairmass = AKTh2.Mag();

                            if (TMath::Abs(125 - AKTjetpairmass) < TMath::Abs(125 - AKTjetpair2Mass)) {
                                AKTjetpair2Mass = AKTjetpairmass;
                                pair2jet1entry = akt1entry2;
                                pair2jet2entry = akt2entry2;

                            }
                        }
                    }
                }
            }

	    */
            
            if (LepFlag == 0) {
                LEPeta = Lep_e_eta -> GetValue(pair1jet1entry);
                LEPphi = Lep_e_phi -> GetValue(pair1jet1entry);
                LEPpt = Lep_e_pt -> GetValue(pair1jet1entry);
                LEPcharge = Lep_e_Charge -> GetValue(pair1jet1entry);
            } else if (LepFlag == 1) {
                LEPeta = Lep_mu_eta -> GetValue(pair1jet1entry);
                LEPphi = Lep_mu_phi -> GetValue(pair1jet1entry);
                LEPpt = Lep_mu_pt -> GetValue(pair1jet1entry);
                LEPcharge = Lep_mu_Charge -> GetValue(pair1jet1entry);
            }
	
	    AKTjet1eta2 = AKTjet_eta -> GetValue(pair1jet2entry);
	    AKTjet1phi2 = AKTjet_phi -> GetValue(pair1jet2entry);
	    AKTjet1pt2 = AKTjet_pt -> GetValue(pair1jet2entry);
	    JetEnergyFix(AKTjet1eta2, AKTjet1pt2, JER);
	   /* 
	    if (entry != 2676 and entry != 3155 and entry != 3515 and entry != 3632) {
		JetEnergyFix(AKTjet1eta2, AKTjet1pt2, JER);
	    } else {
		cout << AKTjet1eta2 << "||" << AKTjet1pt2 << endl;
		cout << pair1jet1entry << "||" << pair1jet2entry << endl;
	    }
	    */
	    AKTjet1mass2 = AKTjet_mass -> GetValue(pair1jet2entry);
	    AKTjet1charge2 = AKTjet_Charge -> GetValue(pair1jet2entry);

	    AKTpairCharge = LEPcharge * AKTjet1charge2;

	    AKTjet2eta1 = AKTjet_eta -> GetValue(pair2jet1entry);
	    AKTjet2phi1 = AKTjet_phi -> GetValue(pair2jet1entry);
	    AKTjet2pt1 = AKTjet_pt -> GetValue(pair2jet1entry);
	    JetEnergyFix(AKTjet2eta1, AKTjet2pt1, JER);
	    AKTjet2mass1 = AKTjet_mass -> GetValue(pair2jet1entry);
        	
	    AKTjet2eta2 = AKTjet_eta -> GetValue(pair2jet2entry);
            AKTjet2phi2 = AKTjet_phi -> GetValue(pair2jet2entry);
            AKTjet2pt2 = AKTjet_pt -> GetValue(pair2jet2entry);
            JetEnergyFix(AKTjet2eta2, AKTjet2pt2, JER);
            AKTjet2mass2 = AKTjet_mass -> GetValue(pair2jet2entry);

            AKTjet1BTag2 = AKTjet_BTag -> GetValue(pair1jet2entry);
            AKTjet2BTag1 = AKTjet_BTag -> GetValue(pair2jet1entry);
            AKTjet2BTag2 = AKTjet_BTag -> GetValue(pair2jet2entry);
	
            AKTjet1TauTag2 = AKTjet_TauTag -> GetValue(pair1jet2entry);
            AKTjet2TauTag1 = AKTjet_TauTag -> GetValue(pair2jet1entry);
            AKTjet2TauTag2 = AKTjet_TauTag -> GetValue(pair2jet2entry);

/*
            bool AKT1jet1flag = false;
            bool AKT1jet2flag = false;
            bool AKT2jet1flag = false;
            bool AKT2jet2flag = false;

            bool AKT1jet1edgeflag = false;
            bool AKT1jet2edgeflag = false;
            bool AKT2jet1edgeflag = false;
            bool AKT2jet2edgeflag = false;
*/
            bool disCutFlag = false;
            bool TauFlag = false;
/*
	    jet1DeltaR = 100;
	    jet2DeltaR = 100;
	    jet3DeltaR = 100;
	    jet4DeltaR = 100;

	    Double_t jet1DeltaPhi;
	    Double_t jet2DeltaPhi;
	    Double_t jet3DeltaPhi;
	    Double_t jet4DeltaPhi;

	    Int_t jet1entry;
	    Int_t jet2entry;
	    Int_t jet3entry;
	    Int_t jet4entry;
*/
	    /*
	    if (AKTjetpair1Mass > 0 and AKTjetpair2Mass >= 0 and AKTjetpair1Mass < 160 and AKTjetpair2Mass < 160 and TMath::Abs(AKTjetpair1Mass-AKTjetpair2Mass) <= 40) {
  	        disCutFlag = true;
	    }
	    */

	    /*
            if (AKTjet1TauTag1 != 0 and AKTjet1TauTag2 != 0) {
                TauFlag = true;
            }
	    */
            Double_t AKTpt[] = {
                LEPpt,
                AKTjet1pt2,
                AKTjet2pt1,
                AKTjet2pt2
            };
	    /*
            Double_t deltaR[] = {
                jet1DeltaR,
                jet2DeltaR,
                jet3DeltaR,
                jet4DeltaR
            };
*/
            Double_t AKTeta[] = {
                LEPeta,
                AKTjet1eta2,
                AKTjet2eta1,
                AKTjet2eta2
            };

            Double_t AKTphi[] = {
                LEPphi,
                AKTjet1phi2,
                AKTjet2phi1,
                AKTjet2phi2
            };

            Double_t AKTmass[] = {
                (1 - LepFlag) * m_e + LepFlag * m_mu,
                AKTjet1mass2,
                AKTjet2mass1,
                AKTjet2mass2
	    	};

            Int_t AKTBTag[] = {
                AKTjet1BTag2,
                AKTjet2BTag1,
                AKTjet2BTag2
            };
/*
	    if (AKTBTag[0] >= 4//AKTBTag[0]%2==1//){
	        BDTjet1BTag1 = 1;
	    } else {
	        BDTjet1BTag1 = 0;
	    }
*/
	    if (AKTBTag[0] >= 4/*AKTBTag[1]%2==1*/){
	        BDTjet1BTag2 = 1;
	    } else {
	        BDTjet1BTag2 = 0;
	    }

	    if (AKTBTag[1] >= 4/*AKTBTag[3]%2==1*/){
	        BDTjet2BTag1 = 1;
	    } else {
	        BDTjet2BTag1 = 0;
	    }

	    if (AKTBTag[2] >= 4/*AKTBTag[3]%2==1*/){
	        BDTjet2BTag2 = 1;
	    } else {
	        BDTjet2BTag2 = 0;
	    }

	    BDTjet1TauTag2 = AKTjet1TauTag2;
	    BDTjet2TauTag1 = AKTjet2TauTag1;
	    BDTjet2TauTag2 = AKTjet2TauTag2;

	    NBTag = /*BDTjet1BTag1 + BDTjet1BTag2 +*/ BDTjet2BTag1 + BDTjet2BTag2;
	    NTauTag = BDTjet1TauTag2 /*+ BDTjet2TauTag1 + BDTjet2TauTag2*/;

	    BDTlep1charge1 = LEPcharge;
	    BDTjet1charge2 = AKTjet1charge2;
	    BDTpairCharge = BDTlep1charge1 * BDTjet1charge2;
/*
	    //calculate angle
            AKTjet1theta1 = 2 *atan(exp(-AKTeta[0]));
            AKTjet1theta2 = 2 *atan(exp(-AKTeta[1]));
            AKTjet2theta1 = 2 *atan(exp(-AKTeta[2]));
            AKTjet2theta2 = 2 *atan(exp(-AKTeta[3]));
*/           
	    //if (AKT1jet1flag == true and AKT1jet2flag == true and AKT2jet1flag == true and AKT2jet2flag == true){
	    AKT1jet1.SetPtEtaPhiM(AKTpt[0], AKTeta[0], AKTphi[0], AKTmass[0]);
	    AKT1jet2.SetPtEtaPhiM(AKTpt[1], AKTeta[1], AKTphi[1], AKTmass[1]);
	    AKT2jet1.SetPtEtaPhiM(AKTpt[2], AKTeta[2], AKTphi[2], AKTmass[2]);
	    AKT2jet2.SetPtEtaPhiM(AKTpt[3], AKTeta[3], AKTphi[3], AKTmass[3]);

	    AKTh1 = AKT1jet1 + AKT1jet2;
	    AKTh2 = AKT2jet1 + AKT2jet2;
	     
	    AKTjetpair1Mass = AKTh1.Mag();
	    AKTjetpair2Mass = AKTh2.Mag();

	    TLorentzVector AKTdiH = AKTh1 + AKTh2;

	    
	    if (AKTh2.Mag() < 50 and AKTh2.Mag() > 0) {
		cout << "beta1 = " << AKTeta[2] << endl;
		cout << "bphi1 = " << AKTphi[2] << endl;
		cout << "bpt1 = " << AKTpt[2] << endl;
		cout << "beta2 = " << AKTeta[3] << endl;
		cout << "bphi2 = " << AKTphi[3] << endl;
		cout << "bpt2 = " << AKTpt[3] << endl;
		cout << "entry = " << entry<< endl;
		/*
		for (Int_t TauIndex = 0; TauIndex < nParticle; TauIndex++) {
		    Int_t TauId = PID -> GetValue(TauIndex);
		    if (TauId == 15 or TauId == -15){
		        Tau1eta = GenParticleEta -> GetValue(TauIndex);
		        Tau1phi = GenParticlePhi -> GetValue(TauIndex);
		        Tau1pt = GenParticlePt -> GetValue(TauIndex);
			    cout << "TauID:::" << "TauEta = " << Tau1eta << "; TauPhi = " << Tau1phi << "; TauPt = " << Tau1pt << endl; 
		    }
		}
		cout << "BBbar mass = " << AKTh2.Mag() << endl;

		cout << "lep1eta1 = " << AKTeta[0] << endl;
		cout << "lep1phi1 = " << AKTphi[0] << endl;
	:	cout << "lep1pt1 = " << AKTpt[0] << endl;
		cout << "jet1eta2 = " << AKTeta[1] << endl;
		cout << "jet1phi2 = " << AKTphi[1] << endl;
		cout << "jet1pt2 = " << AKTpt[1] << endl;
		*/

	    }
	    if (AKTjetpair2Mass < 90 or AKTjetpair2Mass > 130){
	        continue;
	    } 
	    if (NBTag == 2 and NTauTag == 1 and BDTpairCharge == -1 /*and AKTjetpair2Mass > 50*/ and LepFlag != 2  /* and disCutFlag == true and TauFlag == true*/){

	 	AKTjetMass1 -> Fill(AKTjetpair1Mass);
	   	AKTjetMass2 -> Fill(AKTjetpair2Mass);
		AKTjetPairs -> Fill(AKTjetpair1Mass, AKTjetpair2Mass);

		BDTlep1pt1 = AKTpt[0];
		BDTjet1pt2 = AKTpt[1];
		BDTjet2pt1 = AKTpt[2];
		BDTjet2pt2 = AKTpt[3];
		BDTlep1eta1 = AKTeta[0];
		BDTjet1eta2 = AKTeta[1];
		BDTjet2eta1 = AKTeta[2];
		BDTjet2eta2 = AKTeta[3];
		BDTlep1phi1 = AKTphi[0];
		BDTjet1phi2 = AKTphi[1];
		BDTjet2phi1 = AKTphi[2];
		BDTjet2phi2 = AKTphi[3];

		BDTNjets = nAKTjet;
		BDThiggs1pt = AKTh1.Pt();
		BDThiggs1eta = AKTh1.Eta();
		BDThiggs1phi = AKTh1.Phi();
		BDThiggs1invm = AKTh1.Mag();
	
		BDThiggs2pt = AKTh2.Pt();
		BDThiggs2eta = AKTh2.Eta();
		BDThiggs2phi = AKTh2.Phi();
		BDThiggs2invm = AKTh2.Mag();

		BDThiggsAngle = AKTh1.Angle(AKTh2.Vect());

	        BDThiggsdeltaEta = TMath::Abs(BDThiggs1eta - BDThiggs2eta);
	        BDThiggsdeltaPhi = TMath::Abs(BDThiggs1phi - BDThiggs2phi);

	        BDTdihiggspt = AKTdiH.Pt();
	        BDTdihiggseta = AKTdiH.Eta();
	        BDTdihiggsphi = AKTdiH.Phi();
	        BDTdihiggsinvm = AKTdiH.Mag();

  	        DiHiggsMETDeltaTheta = TMath::Abs(2*TMath::ATan(TMath::Exp(-BDTdihiggseta))-2*TMath::ATan(TMath::Exp(-BDTMissingETEta)));

    	        NBTag = /*BDTjet1BTag1 + BDTjet1BTag2 +*/ BDTjet2BTag1 + BDTjet2BTag2;
	        NTauTag = BDTjet1TauTag2 /*+ BDTjet2TauTag1 + BDTjet2TauTag2*/;
 
		//cout << "BDTjet1charge1 = " << BDTjet1charge1 << "BDTjet1charge2 = " << BDTjet1charge2 << endl;
		//cout << "LepFlag = " << LepFlag << "; pair1jet1entry = " << pair1jet1entry << "; pair1jet2entry = " << pair1jet2entry << endl;
		
		if (BDThiggs2invm > 0){
	    	    tree_BDT -> Fill();
		}
	    }
	    
        }
    }	    
    file_sig -> Close();

    //file_BDT -> cd();
    //tree_BDT -> Write();
    //file_BDT -> Close();

    //printf("Time taken: %.2fs\n", (Double_t)(clock() - tStart)/CLOCKS_PER_SEC);
    cout << "Create dataset for Boosting Decision Tress training...." << endl;
    cout << "Done" << endl;
}
void Pairing_tau_mix_6TeV(const char *inputSigFile,
     const char *inputBkg1File,
     	 const char *inputBkg2File,
	     const char *inputBkg3File,
	    	 /*const char *inputBkg4File,*/
	            const char *outputFile,
		    	const char *inputFileForJES,
			    const char *outputFileForJES) {
	
    struct timespec begin, end;

    clock_gettime(CLOCK_REALTIME, &begin);

    //Initiation
    gSystem -> Load("libDelphes.so");
    gROOT -> ProcessLine("gErrorIgnoreLevel = kWarning;");
    gErrorIgnoreLevel = kWarning;

    //gStyle->SetCanvasPreferGL(1);
    //TGaxis::SetMaxDigits(3);
    //Calibration
    Double_t JER[10][10];
    cout << endl << "Calling Calibration algorithm..." << endl;
    Calibration(inputFileForJES, outputFileForJES, JER);

    TFile *output = new TFile(outputFile, "RECREATE");
    TTree *tree_output = new TTree("tree_output", "Delphes");
    TTree *tree_BDT_sig = new TTree("tree_BDT_sig", "sig");
    TTree *tree_BDT_bkg1 = new TTree("tree_BDT_bkg1", "bkg1");
    TTree *tree_BDT_bkg2 = new TTree("tree_BDT_bkg2", "bkg2");
    TTree *tree_BDT_bkg3 = new TTree("tree_BDT_bkg3", "bkg3");
    //TTree *tree_BDT_bkg4 = new TTree("tree_BDT_bkg4", "bkg4");
    
    THStack* JetPair1 = new THStack("JetPair1", ";Tau-jets Pair Invariant Mass [GeV];Events");
    THStack* JetPair2 = new THStack("JetPair2", ";b-jets Pair Invariant Mass [GeV];Events");
    THStack* JetPair2DBDT = new THStack("JetPair2DBDT", ";Tau-jets Pair Invariant Mass [GeV];b-jets Pair Invariant Mass [GeV]");

    TH1D *AKTjetMass1_sig = new TH1D("AKTjetMass1_sig", "Anti_KTjet leading jets pair invariant mass", 16, 0, 160);
    TH1D *AKTjetMass2_sig = new TH1D("AKTjetMass2_sig", "Anti_KTjet sub-leading jets pair invariant mass", 16, 0, 160); 

    TH1D *AKTjetMass1_bkg1 = new TH1D("AKTjetMass1_bkg1", "Anti_KTjet leading jets pair invariant mass", 16, 0, 160);
    TH1D *AKTjetMass2_bkg1 = new TH1D("AKTjetMass2_bkg1", "Anti_KTjet sub-leading jets pair invariant mass", 16, 0, 160); 


    TH1D *AKTjetMass1_bkg2 = new TH1D("AKTjetMass1_bkg2", "Anti_KTjet leading jets pair invariant mass", 16, 0, 160);
    TH1D *AKTjetMass2_bkg2 = new TH1D("AKTjetMass2_bkg2", "Anti_KTjet sub-leading jets pair invariant mass", 16, 0, 160); 


    TH1D *AKTjetMass1_bkg3 = new TH1D("AKTjetMass1_bkg3", "Anti_KTjet leading jets pair invariant mass", 16, 0, 160);
    TH1D *AKTjetMass2_bkg3 = new TH1D("AKTjetMass2_bkg3", "Anti_KTjet sub-leading jets pair invariant mass", 16, 0, 160); 
/*
    TH1D *AKTjetMass1_bkg4 = new TH1D("AKTjetMass1_bkg4", "Anti_KTjet leading jets pair invariant mass", 11, 50, 160);
    TH1D *AKTjetMass2_bkg4 = new TH1D("AKTjetMass2_bkg4", "Anti_KTjet sub-leading jets pair invariant mass", 11, 50, 160); 
*/

    UInt_t nbins = 160;
    UInt_t mmin = 0;
    UInt_t mmax = 160;
    UInt_t bins = 16;

    TH2F *AKTjetPairsMass_sigBDT = new TH2F("AKTjetPairMass_sigBDT", "Anti_KTjet diHiggs invariant mass", nbins, mmin, 160, nbins, mmin, 160);
    TH2F *AKTjetPairsMass_bkg1BDT = new TH2F("AKTjetPairMass_bkg1BDT", "Anti_KTjet diHiggs invariant mass", nbins, mmin, 160, nbins, mmin, 160);
    TH2F *AKTjetPairsMass_bkg2BDT = new TH2F("AKTjetPairMass_bkg2BDT", "Anti_KTjet diHiggs invariant mass", nbins, mmin, 160, nbins, mmin, 160);
    TH2F *AKTjetPairsMass_bkg3BDT = new TH2F("AKTjetPairMass_bkg3BDT", "Anti_KTjet diHiggs invariant mass", nbins, mmin, 160, nbins, mmin, 160);

    Pairing_w_JES(inputSigFile, JER, AKTjetMass1_sig, AKTjetMass2_sig, AKTjetPairsMass_sigBDT, tree_BDT_sig, output); 
    Pairing_w_JES(inputBkg1File, JER, AKTjetMass1_bkg1, AKTjetMass2_bkg1, AKTjetPairsMass_bkg1BDT, tree_BDT_bkg1, output); 
    Pairing_w_JES(inputBkg2File, JER, AKTjetMass1_bkg2, AKTjetMass2_bkg2, AKTjetPairsMass_bkg2BDT, tree_BDT_bkg2, output); 
    Pairing_w_JES(inputBkg3File, JER, AKTjetMass1_bkg3, AKTjetMass2_bkg3, AKTjetPairsMass_bkg3BDT, tree_BDT_bkg3, output); 

    Int_t SigEntries = tree_BDT_sig -> GetEntries();
    Int_t bkg1Entries = tree_BDT_bkg1 -> GetEntries();
    Int_t bkg2Entries = tree_BDT_bkg2 -> GetEntries();
    Int_t bkg3Entries = tree_BDT_bkg3 -> GetEntries();

    Double_t weight1 = 0.003483;
    Double_t weight2 = 0.014122683;
    Double_t weight3 = 0.008601458;
    Double_t weight4 = 0.00734;

    Double_t SigStrength = SigEntries * weight1/TMath::Sqrt(SigEntries * weight1 + bkg1Entries * weight2 + bkg2Entries * weight3 + bkg3Entries * weight4);

    cout << endl << SigStrength << endl;

    AKTjetMass1_bkg1 -> Scale(weight2);
    AKTjetMass1_bkg1 -> SetFillColor(kAzure+7);
    JetPair1 -> Add(AKTjetMass1_bkg1);
    AKTjetMass1_bkg2 -> Scale(weight3);
    AKTjetMass1_bkg2 -> SetFillColor(kAzure+10);
    JetPair1 -> Add(AKTjetMass1_bkg2);
    AKTjetMass1_bkg3 -> Scale(weight4);
    AKTjetMass1_bkg3 -> SetFillColor(kAzure);
    JetPair1 -> Add(AKTjetMass1_bkg3);
    AKTjetMass1_sig -> Scale(weight1);
    AKTjetMass1_sig -> SetFillColor(kPink-1);
    JetPair1 -> Add(AKTjetMass1_sig);
    
    AKTjetMass2_bkg1 -> Scale(weight2);
    AKTjetMass2_bkg1-> SetFillColor(kAzure+7);
    JetPair2 -> Add(AKTjetMass2_bkg1);
    AKTjetMass2_bkg2 -> Scale(weight3);
    AKTjetMass2_bkg2-> SetFillColor(kAzure+10);
    JetPair2 -> Add(AKTjetMass2_bkg2);
    AKTjetMass2_bkg3 -> Scale(weight4);
    AKTjetMass2_bkg3-> SetFillColor(kAzure);
    JetPair2 -> Add(AKTjetMass2_bkg3);
    AKTjetMass2_sig -> Scale(weight1);
    AKTjetMass2_sig -> SetFillColor(kPink-1);
    JetPair2 -> Add(AKTjetMass2_sig);

    AKTjetPairsMass_sigBDT -> Scale(weight1);
    AKTjetPairsMass_sigBDT -> SetMarkerColor(kPink-1);
    AKTjetPairsMass_sigBDT-> SetMarkerStyle(kFullDotLarge);
    JetPair2DBDT -> Add(AKTjetPairsMass_sigBDT);
    AKTjetPairsMass_bkg1BDT -> Scale(weight2);
    AKTjetPairsMass_bkg1BDT-> SetMarkerColor(kAzure+7);
    AKTjetPairsMass_bkg1BDT-> SetMarkerStyle(kFullDotLarge);
    JetPair2DBDT -> Add(AKTjetPairsMass_bkg1BDT);
    AKTjetPairsMass_bkg2BDT -> Scale(weight3);
    AKTjetPairsMass_bkg2BDT-> SetMarkerColor(kAzure+10);
    AKTjetPairsMass_bkg2BDT-> SetMarkerStyle(kFullDotLarge);
    JetPair2DBDT -> Add(AKTjetPairsMass_bkg2BDT);
    AKTjetPairsMass_bkg3BDT -> Scale(weight4);
    AKTjetPairsMass_bkg3BDT-> SetMarkerColor(kAzure);
    AKTjetPairsMass_bkg3BDT-> SetMarkerStyle(kFullDotLarge);
    JetPair2DBDT -> Add(AKTjetPairsMass_bkg3BDT);


    TCanvas *totalcanvas = new TCanvas("totalcanvas", "Canvas", 1400, 1400, 1400, 1400);
    totalcanvas -> SetWindowSize(1204,1228);
    totalcanvas -> SetCanvasSize(1200,1200);

    JetPair1 -> Draw("HIST");
    //JetPair1 -> SetMaximum(180);
    TH1* AKTjetMass1_sig_clone= new TH1D(*AKTjetMass1_sig);
    AKTjetMass1_sig_clone -> SetLineColor(kPink-1);
    AKTjetMass1_sig_clone -> SetLineWidth(4);
    AKTjetMass1_sig_clone -> Draw("SAME");
    TLegend *legend = new TLegend(0.6, 0.65, 0.85, 0.85);
    legend -> AddEntry(AKTjetMass1_sig, "signal HH", "f");
    legend -> AddEntry(AKTjetMass1_bkg1, "ZH(#tau#tau)", "f");
    legend -> AddEntry(AKTjetMass1_bkg2, "Higgs + qq", "f");
    legend -> AddEntry(AKTjetMass1_bkg3, "HZ(#tau#tau)", "f");
    legend -> AddEntry(AKTjetMass1_sig_clone, "signal HH", "l");
    legend -> SetBorderSize(0);
    legend -> Draw();
    TLatex title1(35, 20, "#bf{#font[72]{Muon Collider Simulation (Delphes)}}");
    title1.SetTextSize(0.04); 
    TLatex latexup1(110, 11.5, "#font[52]{#sqrt{s}} #font[42]{=} #font[52]{6 TeV}");
    latexup1.SetTextSize(0.025); 
    TLatex latexdown1(110, 10, "#font[52]{L} #font[42]{=} #font[52]{4 ab^{-1}}");
    latexdown1.SetTextSize(0.025);
    string init1("#font[52]{#frac{S}{#sqrt{S+B}}} #font[42]{=} #font[52]{");
    string add1 = to_string(SigStrength);
    string end1("}");
    //string end1("} #font[52]{(0#leqm_{H}#leq250)}");
    init1 = init1 + add1.substr(0, 5) + end1;
    const char * latex_1 = init1.c_str();
    TLatex latexSigstrength1(110, 8.5, latex_1);
    latexSigstrength1.SetTextSize(0.025);
    latexSigstrength1.Draw("same");
    latexup1.Draw("same");
    latexdown1.Draw("same");
    title1.Draw("same");
    totalcanvas -> SaveAs("JetPairs1_sig+bkg_tau_mix_6TeV.png");

    JetPair2 -> Draw("HIST");
    TH1* AKTjetMass2_sig_clone = new TH1D(*AKTjetMass2_sig);
    AKTjetMass2_sig_clone -> SetLineColor(kPink-1);  
    AKTjetMass2_sig_clone -> SetLineWidth(4);  
    AKTjetMass2_sig_clone -> Draw("SAME");
    TLegend *legend2 = new TLegend(0.15, 0.65, 0.5, 0.85);
    legend2 -> AddEntry(AKTjetMass2_sig, "signal HH", "f");
    legend2 -> AddEntry(AKTjetMass2_bkg1, "ZH(#tau#tau)", "f");
    legend2 -> AddEntry(AKTjetMass2_bkg2, "Higgs + qq", "f");
    legend2 -> AddEntry(AKTjetMass2_bkg3, "HZ(#tau#tau)", "f");
    legend2 -> AddEntry(AKTjetMass2_sig_clone, "signal HH", "l");
    legend2 -> SetBorderSize(0);
    legend2 -> Draw();
    TLatex title2(35, 20.5, "#bf{#font[72]{Muon Collider Simulation (Delphes)}}");
    title2.SetTextSize(0.04); 
    TLatex latexup2(10, 12.5, "#font[52]{#sqrt{s}} #font[42]{=} #font[52]{6 TeV}");
    latexup2.SetTextSize(0.025); 
    TLatex latexdown2(10, 11, "#font[52]{L} #font[42]{=} #font[52]{4 ab^{-1}}");
    latexdown2.SetTextSize(0.025);
    string init2("#font[52]{#frac{S}{#sqrt{S+B}}} #font[42]{=} #font[52]{");
    string add2 = to_string(SigStrength);
    string end2("}");
    //string end2("} #font[52]{(0#leqm_{H}#leq250)}");
    init2 = init2 + add2.substr(0, 5) + end2;
    const char * latex_2 = init2.c_str();
    TLatex latexSigstrength2(10, 9.5, latex_2);
    latexSigstrength2.SetTextSize(0.025);
    latexSigstrength2.Draw("same");
    latexup2.Draw("same");
    latexdown2.Draw("same");
    title2.Draw("same");
    totalcanvas -> SaveAs("JetPairs2_sig+bkg_tau_mix_6TeV.png");

    totalcanvas -> SetCanvasSize(1300,1200);
    totalcanvas -> SetLogy(0);
    JetPair2DBDT -> Draw("SCAT=20");
    TLegend *legend4 = new TLegend(0.65, 0.15, 0.875, 0.35);
    legend4 -> AddEntry(AKTjetPairsMass_sigBDT, "signal HH data", "p");
    legend4 -> AddEntry(AKTjetPairsMass_bkg1BDT, "ZH(#tau#tau)", "p");
    legend4 -> AddEntry(AKTjetPairsMass_bkg2BDT, "Higgs + qq", "p");
    legend4 -> AddEntry(AKTjetPairsMass_bkg3BDT, "HZ(#tau#tau)", "p");
    legend4 -> SetBorderSize(0);
    legend4 -> Draw("same");
    totalcanvas -> cd();
    TLatex title4(40, 161, "#bf{#font[72]{Muon Collider Simulation (Delphes)}}");
    title4.SetTextSize(0.04);
    TLatex latexUp4(120, 155, "#font[52]{#sqrt{s}} #font[42]{=} #font[52]{6 TeV}");
    TLatex latexDown4(120, 145, "#font[52]{L} #font[42]{=} #font[52]{4 ab^{-1}}");
    string init3("#font[52]{#frac{S}{#sqrt{S+B}}} #font[42]{=} #font[52]{");
    string add3 = to_string(SigStrength);
    string end3("}");
    //string end3("} #font[52]{(0#leqm_{H}#leq250)}");
    init3 = init3 + add3.substr(0, 5) + end3;
    const char * latex3 = init3.c_str();
    TLatex latexSigStrength3(120, 135, latex3);
    latexUp4.SetTextSize(0.025);
    latexDown4.SetTextSize(0.025);
    latexSigStrength3.SetTextSize(0.025);
    latexUp4.Draw("same");
    latexDown4.Draw("same");
    latexSigStrength3.Draw("same");
    title4.Draw("same");
    totalcanvas -> SaveAs("JetPair2D_sig+bkgBDT_tau_mix_6TeV.png");
    
    output -> cd();

    JetPair1 -> Write();
    JetPair2 -> Write();

    tree_BDT_sig -> Write();
    tree_BDT_bkg1 -> Write();
    tree_BDT_bkg2 -> Write();
    tree_BDT_bkg3 -> Write();
    //tree_output -> Write();
    output -> Close();
    //train BDT model
    cout << "Initating TMVA for Multi-Variate Analysis...." << endl;
    //apply BDT model
    cout << "Start Boosting Decision Trees application..." << endl;

    clock_gettime(CLOCK_REALTIME, &end);
    long seconds = end.tv_sec - begin.tv_sec;
    long nanoseconds = end.tv_nsec - begin.tv_nsec;
    double elapsed = seconds + nanoseconds*1e-9;
    printf("Time taken: %.3f seconds.\n",elapsed);

    cout << "All done! Exit..." << endl;

}
