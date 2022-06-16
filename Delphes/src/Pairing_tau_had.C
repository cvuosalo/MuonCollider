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

    if (0 < AKTjetTheta and AKTjetTheta <= 0.3) {
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

void CalibrationR02(const char *inputFile,
    const char *outputFile,
        Double_t JER[10][10]) { 
    //Initiation
    gSystem -> Load("libDelphes.so");
    TFile *file_sig = new TFile(inputFile);
    TFile *output = new TFile(outputFile, "recreate");
    TTree *tree_sig = (TTree *) file_sig -> Get("Delphes");
    TTree *tree_output = new TTree("tree_output", "Delphes");

    cout << endl << "Initiating Calibration..." << endl;

    TLeaf *AKTjet_size = tree_sig -> GetLeaf("AKTR02jet_size");
    TLeaf *AKTjet_eta = tree_sig -> GetLeaf("AKTR02jet.Eta");
    TLeaf *AKTjet_phi = tree_sig -> GetLeaf("AKTR02jet.Phi");
    TLeaf *AKTjet_pt = tree_sig -> GetLeaf("AKTR02jet.PT");
    TLeaf *AKTjet_mass = tree_sig -> GetLeaf("AKTR02jet.Mass");

    TLeaf *GenJet_size = tree_sig -> GetLeaf("GenR02Jet_size");
    TLeaf *GenJet_eta = tree_sig -> GetLeaf("GenR02Jet.Eta");
    TLeaf *GenJet_phi = tree_sig -> GetLeaf("GenR02Jet.Phi");
    TLeaf *GenJet_pt = tree_sig -> GetLeaf("GenR02Jet.PT");
    TLeaf *GenJet_mass = tree_sig -> GetLeaf("GenR02Jet.Mass");

    Int_t nEntries = tree_sig -> GetEntries();

    TH1D *AKTGenMass1Comp = new TH1D("AKTGenMass1Comp", "AKTGenMass1Comp", 50, -1, 2);
    TH2D *AKTjetPT_Theta = new TH2D("AKTjetPT_Theta", "AKTjetPT_Theta", 30, 0.1482, 3, 30, 0, 400);
    TH2D *GenJetPT_Theta = new TH2D("GenJetPT_Theta", "GenJetPT_Theta", 30, 0.1482, 3, 30, 0, 400);

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

    TLorentzVector AKTjet;
    TLorentzVector GenJet;

    Double_t MatchedGenJetEta;
    Double_t MatchedGenJetPhi;

    //Double_t JER[10][10];

    Double_t matchJetCnt[10][10];
    for (Int_t ThetaEntry = 0; ThetaEntry < 10; ThetaEntry++) {
        for (Int_t PTentry = 0; PTentry < 10; PTentry++) {
            matchJetCnt[ThetaEntry][PTentry] = 0;
            JER[ThetaEntry][PTentry] = 0;
        }
    }
    bool MuonTagging = false;

    Double_t JERtmp;
    Double_t DeltaPhi;
    Double_t DeltaR;
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

        if (nAKTjet >= 6) {
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

                if (DeltaR < 0.2) {
                    
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

    Int_t nEntries = tree_sig -> GetEntries();

    TH1D *AKTGenMass1Comp = new TH1D("AKTGenMass1Comp", "AKTGenMass1Comp", 50, -1, 2);
    TH2D *AKTjetPT_Theta = new TH2D("AKTjetPT_Theta", "AKTjetPT_Theta", 30, 0.1482, 3, 30, 0, 400);
    TH2D *GenJetPT_Theta = new TH2D("GenJetPT_Theta", "GenJetPT_Theta", 30, 0.1482, 3, 30, 0, 400);

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

    TLorentzVector AKTjet;
    TLorentzVector GenJet;

    Double_t MatchedGenJetEta;
    Double_t MatchedGenJetPhi;

    //Double_t JER[10][10];

    Double_t matchJetCnt[10][10];
    for (Int_t ThetaEntry = 0; ThetaEntry < 10; ThetaEntry++) {
        for (Int_t PTentry = 0; PTentry < 10; PTentry++) {
            matchJetCnt[ThetaEntry][PTentry] = 0;
            JER[ThetaEntry][PTentry] = 0;
        }
    }
    bool MuonTagging = false;

    Double_t JERtmp;
    Double_t DeltaPhi;
    Double_t DeltaR;
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

        if (nAKTjet >= 6) {
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
    	Double_t JERR02[10][10],
	    TH1D* AKTjetMass1,
	        TH1D* AKTjetMass2,
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

    TLeaf *AKTjet_size = tree_sig -> GetLeaf("AKTjet_size");
    TLeaf *AKTjet_eta = tree_sig -> GetLeaf("AKTjet.Eta");
    TLeaf *AKTjet_phi = tree_sig -> GetLeaf("AKTjet.Phi");
    TLeaf *AKTjet_pt = tree_sig -> GetLeaf("AKTjet.PT");
    TLeaf *AKTjet_mass = tree_sig -> GetLeaf("AKTjet.Mass");
    TLeaf *AKTjet_BTag = tree_sig -> GetLeaf("AKTjet.BTag");
    TLeaf *AKTjet_TauTag = tree_sig -> GetLeaf("AKTjet.TauTag");
    TLeaf *AKTjet_Charge = tree_sig -> GetLeaf("AKTjet.Charge");
/*
    TLeaf *AKTR02jet_size = tree_sig -> GetLeaf("AKTR01jet_size");
    TLeaf *AKTR02jet_eta = tree_sig -> GetLeaf("AKTR01jet.Eta");
    TLeaf *AKTR02jet_phi = tree_sig -> GetLeaf("AKTR01jet.Phi");
    TLeaf *AKTR02jet_pt = tree_sig -> GetLeaf("AKTR01jet.PT");
    TLeaf *AKTR02jet_mass = tree_sig -> GetLeaf("AKTR01jet.Mass");
    TLeaf *AKTR02jet_BTag = tree_sig -> GetLeaf("AKTR01jet.BTag");
    TLeaf *AKTR02jet_TauTag = tree_sig -> GetLeaf("AKTR01jet.TauTag");
    TLeaf *AKTR02jet_Charge = tree_sig -> GetLeaf("AKTR01jet.Charge");
*/
    TLeaf *AKTR02jet_size = tree_sig -> GetLeaf("AKTR02jet_size");
    TLeaf *AKTR02jet_eta = tree_sig -> GetLeaf("AKTR02jet.Eta");
    TLeaf *AKTR02jet_phi = tree_sig -> GetLeaf("AKTR02jet.Phi");
    TLeaf *AKTR02jet_pt = tree_sig -> GetLeaf("AKTR02jet.PT");
    TLeaf *AKTR02jet_mass = tree_sig -> GetLeaf("AKTR02jet.Mass");
    TLeaf *AKTR02jet_BTag = tree_sig -> GetLeaf("AKTR02jet.BTag");
    TLeaf *AKTR02jet_TauTag = tree_sig -> GetLeaf("AKTR02jet.TauTag");
    TLeaf *AKTR02jet_Charge = tree_sig -> GetLeaf("AKTR02jet.Charge");

    TLeaf *MissingETMET = tree_sig -> GetLeaf("MissingET.MET");
    TLeaf *MissingETEta = tree_sig -> GetLeaf("MissingET.Eta");
    TLeaf *MissingETPhi = tree_sig -> GetLeaf("MissingET.Phi");

    Int_t nEntries = tree_sig -> GetEntries();

    Int_t hasJet = 0;
    Int_t TauSelection = 0;
    Int_t B2 = 0;
    Int_t BSelection = 0;
    Double_t AKTjet1eta1;
    Double_t AKTjet1theta1;
    Double_t AKTjet1phi1;
    Double_t AKTjet1pt1;
    Double_t AKTjet1mass1;
    Int_t AKTjet1charge1;
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

    Int_t AKTjetpairCharge;

    Int_t AKTjet1BTag1;
    Int_t AKTjet1BTag2;
    Int_t AKTjet2BTag1;
    Int_t AKTjet2BTag2;

    Int_t AKTjet1TauTag1;
    Int_t AKTjet1TauTag2;
    Int_t AKTjet2TauTag1;
    Int_t AKTjet2TauTag2;

    Int_t pair1jet1entry;
    Int_t pair1jet2entry;
    Int_t pair2jet1entry;
    Int_t pair2jet2entry;
/*
    Int_t AKTjet2entry1;
    Int_t AKTjet2entry2;
*/
    Double_t diHiggsdis;

    Double_t AKTjetpair1Mass;
    Double_t AKTjetpair2Mass;

    Int_t AKTjetpaircnt;

    TLorentzVector AKTh1;
    TLorentzVector AKTh2;
    TLorentzVector AKT1jet1;
    TLorentzVector AKT1jet2;
    TLorentzVector AKT2jet1;
    TLorentzVector AKT2jet2;

    //BDT variable
    Int_t BDTNjets;

    Double_t BDTjet1pt1;
    Double_t BDTjet1pt2;
    Double_t BDTjet2pt1;
    Double_t BDTjet2pt2;

    Double_t BDTjet1eta1;
    Double_t BDTjet1eta2;
    Double_t BDTjet2eta1;
    Double_t BDTjet2eta2;

    Double_t BDTjet1phi1;
    Double_t BDTjet1phi2;
    Double_t BDTjet2phi1;
    Double_t BDTjet2phi2;

    Int_t BDTjet1charge1;
    Int_t BDTjet1charge2;

    Int_t BDTjetpairCharge;

    Int_t BDTjet1BTag1;
    Int_t BDTjet1BTag2;
    Int_t BDTjet2BTag1;
    Int_t BDTjet2BTag2;

    Int_t BDTjet1TauTag1;
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
   
   
    Double_t BDTDiTauMETDeltaR;
    Double_t dphi;
    Double_t DiHiggsMETDeltaTheta;

    Int_t NBTag;
    Int_t NTauTag;

    tree_BDT -> Branch("BDTNjets", &BDTNjets, "BDTNjets/I");
    tree_BDT -> Branch("NBTag", &NBTag, "NBTag/I");
    tree_BDT -> Branch("NTauTag", &NTauTag, "NTauTag/I");

    tree_BDT -> Branch("BDTjet1BTag1", &BDTjet1BTag1, "BDTjet1BTag1/I");
    tree_BDT -> Branch("BDTjet1BTag2", &BDTjet1BTag2, "BDTjet1BTag2/I");
    tree_BDT -> Branch("BDTjet2BTag1", &BDTjet2BTag1, "BDTjet2BTag1/I");
    tree_BDT -> Branch("BDTjet2BTag2", &BDTjet2BTag2, "BDTjet2BTag2/I");

    tree_BDT -> Branch("BDTjet1TauTag1", &BDTjet1TauTag1, "BDTjet1TauTag1/I");
    tree_BDT -> Branch("BDTjet1TauTag2", &BDTjet1TauTag2, "BDTjet1TauTag2/I");
    tree_BDT -> Branch("BDTjet2TauTag1", &BDTjet2TauTag1, "BDTjet2TauTag1/I");
    tree_BDT -> Branch("BDTjet2TauTag2", &BDTjet2TauTag2, "BDTjet2TauTag2/I");

    tree_BDT -> Branch("BDTjet1pt1", &BDTjet1pt1, "BDTjet1pt1/D");
    tree_BDT -> Branch("BDTjet1pt2", &BDTjet1pt2, "BDTjet1pt2/D");
    tree_BDT -> Branch("BDTjet2pt1", &BDTjet2pt1, "BDTjet2pt1/D");
    tree_BDT -> Branch("BDTjet2pt2", &BDTjet2pt2, "BDTjet2pt2/D");

    tree_BDT -> Branch("BDTjet1eta1", &BDTjet1eta1, "BDTjet1eta1/D");
    tree_BDT -> Branch("BDTjet1eta2", &BDTjet1eta2, "BDTjet1eta2/D");
    tree_BDT -> Branch("BDTjet2eta1", &BDTjet2eta1, "BDTjet2eta1/D");
    tree_BDT -> Branch("BDTjet2eta2", &BDTjet2eta2, "BDTjet2eta2/D");

    tree_BDT -> Branch("BDTjet1phi1", &BDTjet1phi1, "BDTjet1phi1/D");
    tree_BDT -> Branch("BDTjet1phi2", &BDTjet1phi2, "BDTjet1phi2/D");
    tree_BDT -> Branch("BDTjet2phi1", &BDTjet2phi1, "BDTjet2phi1/D");
    tree_BDT -> Branch("BDTjet2phi2", &BDTjet2phi2, "BDTjet2phi2/D");

    tree_BDT -> Branch("BDTjet1charge1", &BDTjet1charge1, "BDTjet1charge1/I");
    tree_BDT -> Branch("BDTjet1charge2", &BDTjet1charge2, "BDTjet1charge2/I");

    tree_BDT -> Branch("BDTjetpairCharge", &BDTjetpairCharge, "BDTjetpairCharge/I");

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
    
    tree_BDT -> Branch("BDTDiTauMETDeltaR", &BDTDiTauMETDeltaR, "BDTDiTauMETDeltaR/D");
    tree_BDT -> Branch("DiHiggsMETDeltaTheta", &DiHiggsMETDeltaTheta, "DiHiggsMETDeltaTheta/D");
    tree_BDT -> Branch("dphi", &dphi, "dphi/D");

    tree_BDT -> Branch("BDTMissingETMET", &BDTMissingETMET, "BDTMissingETMET/D");
    tree_BDT -> Branch("BDTMissingETEta", &BDTMissingETEta, "BDTMissingETEta/D");
    tree_BDT -> Branch("BDTMissingETPhi", &BDTMissingETPhi, "BDTMissingETPhi/D");

    tree_BDT -> Branch("BDTdihiggspt", &BDTdihiggspt, "BDTdihiggspt/D");
    tree_BDT -> Branch("BDTdihiggseta", &BDTdihiggseta, "BDTdihiggseta/D");
    tree_BDT -> Branch("BDTdihiggsphi", &BDTdihiggsphi, "BDTdihiggsphi/D");
    tree_BDT -> Branch("BDTdihiggsinvm", &BDTdihiggsinvm, "BDTdihiggsinvm/D");

    //run over event entries
    for (Long64_t entry = 0; entry < nEntries; entry++) {
        //Initiation
        tree_sig -> GetEntry(entry);
        //tree_output -> GetEntry(entry);
        AKTjet_size -> GetBranch() -> GetEntry(entry);
        AKTR02jet_size -> GetBranch() -> GetEntry(entry);

        Int_t nAKTjet = AKTjet_size -> GetValue();
        Int_t nAKTR02jet = AKTR02jet_size -> GetValue();
	
        AKTjet_eta -> GetBranch() -> GetEntry(entry);
        AKTjet_phi -> GetBranch() -> GetEntry(entry);
        AKTjet_pt -> GetBranch() -> GetEntry(entry);
        AKTjet_mass -> GetBranch() -> GetEntry(entry);
        AKTjet_BTag -> GetBranch() -> GetEntry(entry);
	AKTjet_TauTag -> GetBranch() -> GetEntry(entry);
	AKTjet_Charge -> GetBranch() -> GetEntry(entry);

        AKTR02jet_eta -> GetBranch() -> GetEntry(entry);
        AKTR02jet_phi -> GetBranch() -> GetEntry(entry);
        AKTR02jet_pt -> GetBranch() -> GetEntry(entry);
        AKTR02jet_mass -> GetBranch() -> GetEntry(entry);
        AKTR02jet_BTag -> GetBranch() -> GetEntry(entry);
	AKTR02jet_TauTag -> GetBranch() -> GetEntry(entry);
	AKTR02jet_Charge -> GetBranch() -> GetEntry(entry);

	MissingETMET -> GetBranch() -> GetEntry(entry);
	MissingETEta -> GetBranch() -> GetEntry(entry);
	MissingETPhi -> GetBranch() -> GetEntry(entry);

	BDTMissingETMET = MissingETMET -> GetValue();
        BDTMissingETEta = MissingETEta -> GetValue();
        BDTMissingETPhi = MissingETPhi -> GetValue();

        AKTjetpair1Mass = 0;
        AKTjetpair2Mass = 0;

        pair1jet1entry = 0;
        pair1jet2entry = 0;
        pair2jet1entry = 0;
        pair2jet2entry = 0;

        diHiggsdis = 1000000000;

        AKTjet1BTag1 = 0;
        AKTjet1BTag2 = 0;
        AKTjet2BTag1 = 0;
        AKTjet2BTag2 = 0;
        
        AKTjet1TauTag1 = 0;
        AKTjet1TauTag2 = 0;
        AKTjet2TauTag1 = 0;
        AKTjet2TauTag2 = 0;

        BDTjet1BTag1 = 0;
        BDTjet1BTag2 = 0;
        BDTjet2BTag1 = 0;
        BDTjet2BTag2 = 0;
        
        BDTjet1TauTag1 = 0;
        BDTjet1TauTag2 = 0;
        BDTjet2TauTag1 = 0;
        BDTjet2TauTag2 = 0;

	TVector3 MET;
        TVector3 di_tau;
	Double_t nu_pt=0; 
	Double_t AKTh1Pt=0;
	Double_t Ratio_col=0;
	TLorentzVector Nu;

	bool TauSelectFlag;
        //Pairing up Anti-kt jets
        Double_t AKTjetpair[nAKTR02jet * (nAKTR02jet - 1) / 2][3];
	for (Int_t index = 0; index < (nAKTR02jet*(nAKTR02jet - 1) / 2); index++) {
	    AKTjetpair[index][0] = 0;
	    AKTjetpair[index][1] = 0;
	    AKTjetpair[index][2] = 0;
	}
        AKTjetpaircnt = 0;

        if (nAKTjet >= 3 and nAKTR02jet >= 3 ) {
	    hasJet++;
/*
            for (Int_t aktentry = 0; aktentry < nAKTjet; aktentry++) {
                AKTjetEta = AKTjet_eta->GetValue(aktentry);
                AKTjetPt = AKTjet_pt->GetValue(aktentry);
                AKTjetTheta = 2 * atan(exp(-AKTjetEta));
                AKTjetPT_Theta->Fill(AKTjetTheta, AKTjetPt);
                JetEnergyFix(AKTjetEta, AKTjetPt, JER);
            }
  */  
            //Pairing up leading Tau_Tagged anti-KT jet pair
	    
	    TauSelectFlag = false;
            for (Int_t akt1entry = 0; akt1entry < nAKTR02jet; akt1entry++) {
                for (Int_t akt2entry = akt1entry + 1; akt2entry < nAKTR02jet; akt2entry++) {
                    AKTjet1BTag1 = AKTR02jet_BTag -> GetValue(akt1entry);
                    AKTjet1BTag2 = AKTR02jet_BTag -> GetValue(akt2entry);
                    AKTjet1TauTag1 = AKTR02jet_TauTag -> GetValue(akt1entry);
                    AKTjet1TauTag2 = AKTR02jet_TauTag -> GetValue(akt2entry);
		    AKTjet1charge1 = AKTR02jet_Charge -> GetValue(akt1entry);
		    AKTjet1charge2 = AKTR02jet_Charge -> GetValue(akt2entry);
		    AKTjetpairCharge = AKTjet1charge1 * AKTjet1charge2;

		    AKTjet1eta1 = AKTR02jet_eta -> GetValue(akt1entry);
		    AKTjet1phi1 = AKTR02jet_phi -> GetValue(akt1entry);
		    AKTjet1eta2 = AKTR02jet_eta -> GetValue(akt2entry);
		    AKTjet1phi2 = AKTR02jet_phi -> GetValue(akt2entry);

		    TVector3 TauJet1;
		    TVector3 TauJet2;
		    TauJet1.SetPtEtaPhi(1.0, AKTjet1eta1, AKTjet1phi1);
		    TauJet2.SetPtEtaPhi(1.0, AKTjet1eta2, AKTjet1phi2);
                    if (AKTjet1TauTag1 == 1 and AKTjet1TauTag2 == 1 and AKTjetpairCharge == -1 and AKTjet1BTag1 < 4 and AKTjet1BTag2 < 4/* and TauJet1.DrEtaPhi(TauJet2) <= 3*/ ) {

			//cout << "akt1entry = " << akt1entry << "; akt2entry = " << akt2entry<< endl ;
                        AKTjet1eta1 = AKTR02jet_eta -> GetValue(akt1entry);
                        AKTjet1phi1 = AKTR02jet_phi -> GetValue(akt1entry);
                        AKTjet1pt1 = AKTR02jet_pt -> GetValue(akt1entry);
                        JetEnergyFix(AKTjet1eta1, AKTjet1pt1, JERR02);
                        AKTjet1mass1 = AKTR02jet_mass -> GetValue(akt1entry);
                        AKTjet1charge1 = AKTR02jet_Charge -> GetValue(akt1entry);

                        AKTjet1eta2 = AKTR02jet_eta -> GetValue(akt2entry);
                        AKTjet1phi2 = AKTR02jet_phi -> GetValue(akt2entry);
                        AKTjet1pt2 = AKTR02jet_pt -> GetValue(akt2entry);
                        JetEnergyFix(AKTjet1eta2, AKTjet1pt2, JERR02);
                        AKTjet1mass2 = AKTR02jet_mass -> GetValue(akt2entry);
                        AKTjet1charge2 = AKTR02jet_Charge -> GetValue(akt2entry);

                        AKT1jet1.SetPtEtaPhiM(AKTjet1pt1, AKTjet1eta1, AKTjet1phi1, AKTjet1mass1);
                        AKT1jet2.SetPtEtaPhiM(AKTjet1pt2, AKTjet1eta2, AKTjet1phi2, AKTjet1mass2);
                        AKTh1 = AKT1jet1 + AKT1jet2;
	    
			MET.SetPtEtaPhi(BDTMissingETMET, BDTMissingETEta, BDTMissingETPhi);
			di_tau.SetPtEtaPhi(1.0, AKTh1.Eta(), AKTh1.Phi());
			nu_pt = MET.Dot(di_tau); 
			AKTh1Pt = AKTh1.Pt();
			Ratio_col = TMath::Sqrt((nu_pt + AKTh1Pt)/AKTh1Pt);
			Nu.SetPtEtaPhiM(nu_pt, AKTh1.Eta(), AKTh1.Phi(), 0);
			AKTh1 = AKT1jet1 + AKT1jet2 + Nu;
				
                        AKTjetpairmass = /*Ratio_col **/ AKTh1.Mag();

                        AKTjetpair[AKTjetpaircnt][0] = AKTjetpairmass;
                        AKTjetpair[AKTjetpaircnt][1] = akt1entry;
                        AKTjetpair[AKTjetpaircnt][2] = akt2entry;
                        AKTjetpaircnt++;

                    }
                }
            }

            AKTjetpair1Mass = 0;
            for (Int_t pairindex = 0; pairindex < AKTjetpaircnt /*(nAKTR02jet - 1)*nAKTR02jet/2*/; pairindex++) {
                if (abs(AKTjetpair[pairindex][0] - 125) <= abs(AKTjetpair1Mass - 125)) {
		    TauSelectFlag = true;
                    AKTjetpair1Mass = AKTjetpair[pairindex][0];
                    pair1jet1entry = AKTjetpair[pairindex][1];
                    pair1jet2entry = AKTjetpair[pairindex][2];
            	}
	    }
	    
	    if (pair1jet1entry == 0 and pair1jet2entry == 0) {
		continue;
	    } else {
		TauSelection++;
	    }
	    
	    //cout << "pair1jet1entry = " << pair1jet1entry << "; pair1jet2entry = " << pair1jet2entry ;
	    Int_t BtagSum = 0;
	    Int_t BtagInd = 0;
	    for (Int_t aktentry = 0; aktentry < nAKTjet; aktentry++){
		BtagInd = AKTjet_BTag -> GetValue(aktentry);
		if (BtagInd >=4 ){
		 BtagSum++;
		}
	    }
	    if (BtagSum >= 2){
	        B2++;
	    }
	    // Pairing up the bbbar pair
            AKTjetpair2Mass = 0;
            for (Int_t akt1entry2 = 0; akt1entry2 < nAKTjet; akt1entry2++) {
                AKTjet2eta1 = AKTjet_eta -> GetValue(akt1entry2);
                AKTjet2phi1 = AKTjet_phi -> GetValue(akt1entry2);
                AKTjet2BTag1 = AKTjet_BTag -> GetValue(akt1entry2);
                AKTjet2TauTag1 = AKTjet_TauTag -> GetValue(akt1entry2);

                AKTjet1eta1 = AKTR02jet_eta -> GetValue(pair1jet1entry);
                AKTjet1eta2 = AKTR02jet_eta -> GetValue(pair1jet2entry);
                AKTjet1phi1 = AKTR02jet_phi -> GetValue(pair1jet1entry);
                AKTjet1phi2 = AKTR02jet_phi -> GetValue(pair1jet2entry);

		TVector3 DR1jet1;
		TVector3 DR1jet2;
		TVector3 DR2jet1;
		TVector3 DR2jet2;

		DR1jet1.SetPtEtaPhi(1.0, AKTjet1eta1, AKTjet1phi1);
		DR1jet2.SetPtEtaPhi(1.0, AKTjet1eta2, AKTjet1phi2);
		DR2jet1.SetPtEtaPhi(1.0, AKTjet2eta1, AKTjet2phi1);
                if (/*DR1jet1.DrEtaPhi(DR2jet1) >= 0.4 and DR1jet2.DrEtaPhi(DR2jet1) >= 0.4*/AKTjet2TauTag1 == 0 and AKTjet2BTag1 >= 4) {
                    for (Int_t akt2entry2 = 0; akt2entry2 < nAKTjet; akt2entry2++) {
                	AKTjet2eta2 = AKTjet_eta -> GetValue(akt2entry2);
                	AKTjet2phi2 = AKTjet_phi -> GetValue(akt2entry2);
                	AKTjet2BTag2 = AKTjet_BTag -> GetValue(akt2entry2);
                	AKTjet2TauTag2 = AKTjet_TauTag -> GetValue(akt2entry2);

			DR2jet2.SetPtEtaPhi(1.0, AKTjet2eta2, AKTjet2phi2);
                        if (/*DR1jet1.DrEtaPhi(DR2jet2) >= 0.4 and DR1jet2.DrEtaPhi(DR2jet2) >= 0.4*/ AKTjet2TauTag2 == 0 and akt2entry2 != akt1entry2 and AKTjet2BTag2 >= 4) {

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


                            if (abs(125 - AKTjetpairmass) < abs(125 - AKTjetpair2Mass)) {
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
	    } else {
		BSelection++;
	    }
	    
	    //cout << "; pair2jet1entry = " << pair2jet1entry << "; pair2jet2entry = " << pair2jet2entry << endl;
            AKTjet1eta1 = AKTR02jet_eta -> GetValue(pair1jet1entry);
            AKTjet1phi1 = AKTR02jet_phi -> GetValue(pair1jet1entry);
            AKTjet1pt1 = AKTR02jet_pt -> GetValue(pair1jet1entry);
            JetEnergyFix(AKTjet1eta1, AKTjet1pt1, JERR02);
	    AKTjet1mass1 = AKTR02jet_mass -> GetValue(pair1jet1entry);
	    AKTjet1charge1 = AKTR02jet_Charge -> GetValue(pair1jet1entry);
            
	    AKTjet1eta2 = AKTR02jet_eta -> GetValue(pair1jet2entry);
	    AKTjet1phi2 = AKTR02jet_phi -> GetValue(pair1jet2entry);
	    AKTjet1pt2 = AKTR02jet_pt -> GetValue(pair1jet2entry);
	    JetEnergyFix(AKTjet1eta2, AKTjet1pt2, JERR02);
	    AKTjet1mass2 = AKTR02jet_mass -> GetValue(pair1jet2entry);
	    AKTjet1charge2 = AKTR02jet_Charge -> GetValue(pair1jet2entry);

	    AKTjetpairCharge = AKTjet1charge1 * AKTjet1charge2;

	    AKTjet2eta1 = AKTjet_eta -> GetValue(pair2jet1entry);
	    AKTjet2phi1 = AKTjet_phi -> GetValue(pair2jet1entry);
	    AKTjet2pt1 = AKTjet_pt -> GetValue(pair2jet1entry);
	    JetEnergyFix(AKTjet2eta1, AKTjet2pt1, JER);
	    AKTjet2mass1 = AKTjet_mass -> GetValue(pair2jet1entry);
            
	    AKTjet2eta2 = AKTjet_eta -> GetValue(pair2jet2entry);
            AKTjet2phi2 = AKTjet_phi -> GetValue(pair2jet2entry);
            AKTjet2pt2 = AKTjet_pt -> GetValue(pair2jet2entry);
            AKTjet2mass2 = AKTjet_mass -> GetValue(pair2jet2entry);

	    AKTjet1BTag1 = AKTR02jet_BTag -> GetValue(pair1jet1entry);
            AKTjet1BTag2 = AKTR02jet_BTag -> GetValue(pair1jet2entry);
            AKTjet2BTag1 = AKTjet_BTag -> GetValue(pair2jet1entry);
            AKTjet2BTag2 = AKTjet_BTag -> GetValue(pair2jet2entry);

            AKTjet1TauTag1 = AKTR02jet_TauTag -> GetValue(pair1jet1entry);
            AKTjet1TauTag2 = AKTR02jet_TauTag -> GetValue(pair1jet2entry);
            AKTjet2TauTag1 = AKTjet_TauTag -> GetValue(pair2jet1entry);
            AKTjet2TauTag2 = AKTjet_TauTag -> GetValue(pair2jet2entry);

            bool disCutFlag = false;
            
	    if (AKTjetpair1Mass > 50 and AKTjetpair2Mass >= 0 and AKTjetpair1Mass < 160 and AKTjetpair2Mass < 160 and TMath::Abs(AKTjetpair1Mass-AKTjetpair2Mass) <= 40) {
  	        disCutFlag = true;
	    }

            Double_t AKTpt[] = {
                AKTjet1pt1,
                AKTjet1pt2,
                AKTjet2pt1,
                AKTjet2pt2
            };

            Double_t AKTeta[] = {
                AKTjet1eta1,
                AKTjet1eta2,
                AKTjet2eta1,
                AKTjet2eta2
            };

            Double_t AKTphi[] = {
                AKTjet1phi1,
                AKTjet1phi2,
                AKTjet2phi1,
                AKTjet2phi2
            };

            Double_t AKTmass[] = {
                AKTjet1mass1,
                AKTjet1mass2,
                AKTjet2mass1,
                AKTjet2mass2
	    };

            Int_t AKTBTag[] = {
               	AKTjet1BTag1,
                AKTjet1BTag2,
                AKTjet2BTag1,
                AKTjet2BTag2
            };

	    if (AKTBTag[0] >= 4/*AKTBTag[0]%2==1*/){
	        BDTjet1BTag1 = 1;
	    } else {
	        BDTjet1BTag1 = 0;
	    }

	    if (AKTBTag[1] >= 4/*AKTBTag[1]%2==1*/){
	        BDTjet1BTag2 = 1;
	    } else {
	        BDTjet1BTag2 = 0;
	    }

	    if (AKTBTag[2] >= 4/*AKTBTag[3]%2==1*/){
	        BDTjet2BTag1 = 1;
	    } else {
	        BDTjet2BTag1 = 0;
	    }

	    if (AKTBTag[3] >= 4/*AKTBTag[3]%2==1*/){
	        BDTjet2BTag2 = 1;
	    } else {
	        BDTjet2BTag2 = 0;
	    }

	    BDTjet1TauTag1 = AKTjet1TauTag1;
	    BDTjet1TauTag2 = AKTjet1TauTag2;
	    BDTjet2TauTag1 = AKTjet2TauTag1;
	    BDTjet2TauTag2 = AKTjet2TauTag2;

    	    NBTag = /*BDTjet1BTag1 + BDTjet1BTag2 +*/ BDTjet2BTag1 + BDTjet2BTag2;
	    NTauTag = BDTjet1TauTag1 + BDTjet1TauTag2 /*+ BDTjet2TauTag1 + BDTjet2TauTag2*/;

	    BDTjet1charge1 = AKTjet1charge1;
	    BDTjet1charge2 = AKTjet1charge2;
	    BDTjetpairCharge = AKTjetpairCharge;
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

	    MET.SetPtEtaPhi(BDTMissingETMET, BDTMissingETEta, BDTMissingETPhi);
	    di_tau.SetPtEtaPhi(1.0, AKTh1.Eta(), AKTh1.Phi());
	    nu_pt = MET.Dot(di_tau); 
	    AKTh1Pt = AKTh1.Pt();
	    Ratio_col = TMath::Sqrt((nu_pt + AKTh1Pt)/AKTh1Pt);
	    Nu.SetPtEtaPhiM(nu_pt, AKTh1.Eta(), AKTh1.Phi(), 0);
	    AKTh1 = AKT1jet1 + AKT1jet2 + Nu;
	    AKTjetpair1Mass = /*Ratio_col **/ AKTh1.Mag();
	    AKTjetpair2Mass = AKTh2.Mag();

	    TLorentzVector AKTdiH = Ratio_col * AKTh1 + AKTh2;
	    if (NBTag == 2 and NTauTag == 2 and BDTjetpairCharge == -1 /*and disCutFlag == true */){

	 	AKTjetMass1 -> Fill(AKTjetpair1Mass);
	   	AKTjetMass2 -> Fill(AKTjetpair2Mass);
		
		BDTjet1pt1 = AKTpt[0];
		BDTjet1pt2 = AKTpt[1];
		BDTjet2pt1 = AKTpt[2];
		BDTjet2pt2 = AKTpt[3];
		BDTjet1eta1 = AKTeta[0];
		BDTjet1eta2 = AKTeta[1];
		BDTjet2eta1 = AKTeta[2];
		BDTjet2eta2 = AKTeta[3];
		BDTjet1phi1 = AKTphi[0];
		BDTjet1phi2 = AKTphi[1];
		BDTjet2phi1 = AKTphi[2];
		BDTjet2phi2 = AKTphi[3];

		BDTNjets = nAKTjet;
		BDThiggs1pt = /*Ratio_col **/ AKTh1.Pt();
		BDThiggs1eta = AKTh1.Eta();
		BDThiggs1phi = AKTh1.Phi();
		BDThiggs1invm = /*Ratio_col **/ AKTh1.Mag()/**125/75*/;
	
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

		dphi = TVector2::Phi_mpi_pi(BDThiggs1phi - BDTMissingETPhi);
		BDTDiTauMETDeltaR = TMath::Sqrt((BDThiggs1eta - BDTMissingETEta) * (BDThiggs1eta - BDTMissingETEta) + dphi * dphi);
  	        DiHiggsMETDeltaTheta = TMath::Abs(2*TMath::ATan(TMath::Exp(-BDTdihiggseta))-2*TMath::ATan(TMath::Exp(-BDTMissingETEta)));

    	        NBTag = /*BDTjet1BTag1 + BDTjet1BTag2 +*/ BDTjet2BTag1 + BDTjet2BTag2;
	        NTauTag = BDTjet1TauTag1 + BDTjet1TauTag2 /*+ BDTjet2TauTag1 + BDTjet2TauTag2*/;
		//cout << "BDTjetpairCharge = " << BDTjetpairCharge << endl;
	    	tree_BDT -> Fill();
	    }
	    
        }
    }	    
    file_sig -> Close();
    //file_BDT -> cd();
    //tree_BDT -> Write();
    //file_BDT -> Close();
    cout << "hasJet = " << hasJet << endl;
    cout << "TauSelection = " << TauSelection << endl;
    cout << "B2 = " << B2 << endl;
    cout << "BSelection = " << BSelection << endl;
    //printf("Time taken: %.2fs\n", (Double_t)(clock() - tStart)/CLOCKS_PER_SEC);
    cout << "Create dataset for Boosting Decision Tress training...." << endl;
    cout << "Done" << endl;
}
void Pairing_tau_had(const char *inputSigFile,
    const char *inputBkg1File,
    	const char *inputBkg2File,
	    /*const char *inputBkg3File,
	    	const char *inputBkg4File,*/
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
    Double_t JERR02[10][10];
    cout << endl << "Calling Calibration algorithm..." << endl;
    Calibration(inputFileForJES, outputFileForJES, JER);
    CalibrationR02(inputFileForJES, outputFileForJES, JERR02);

    TFile *output = new TFile(outputFile, "RECREATE");
    //TTree *tree_output = new TTree("tree_output", "Delphes");
    TTree *tree_BDT_sig = new TTree("tree_BDT_sig", "sig");
    TTree *tree_BDT_bkg1 = new TTree("tree_BDT_bkg1", "bkg1");
    TTree *tree_BDT_bkg2 = new TTree("tree_BDT_bkg2", "bkg2");
    //TTree *tree_BDT_bkg3 = new TTree("tree_BDT_bkg3", "bkg3");
    //TTree *tree_BDT_bkg4 = new TTree("tree_BDT_bkg4", "bkg4");
    
    //THStack* JetPair1 = new THStack("JetPair1", ";Leading b-jets Pair Invariant Mass [GeV];Events");
    //THStack* JetPair2 = new THStack("JetPair2", ";Sub-Leading b-jets Pair Invariant Mass [GeV];Events");

    TH1D *AKTjetMass1_sig = new TH1D("AKTjetMass1_sig", "Anti_KTjet leading jets pair invariant mass", 11, 50, 160);
    TH1D *AKTjetMass2_sig = new TH1D("AKTjetMass2_sig", "Anti_KTjet sub-leading jets pair invariant mass", 11, 50, 160); 

    TH1D *AKTjetMass1_bkg1 = new TH1D("AKTjetMass1_bkg1", "Anti_KTjet leading jets pair invariant mass", 11, 50, 160);
    TH1D *AKTjetMass2_bkg1 = new TH1D("AKTjetMass2_bkg1", "Anti_KTjet sub-leading jets pair invariant mass", 11, 50, 160); 


    TH1D *AKTjetMass1_bkg2 = new TH1D("AKTjetMass1_bkg2", "Anti_KTjet leading jets pair invariant mass", 11, 50, 160);
    TH1D *AKTjetMass2_bkg2 = new TH1D("AKTjetMass2_bkg2", "Anti_KTjet sub-leading jets pair invariant mass", 11, 50, 160); 

/*
    TH1D *AKTjetMass1_bkg3 = new TH1D("AKTjetMass1_bkg3", "Anti_KTjet leading jets pair invariant mass", 11, 50, 160);
    TH1D *AKTjetMass2_bkg3 = new TH1D("AKTjetMass2_bkg3", "Anti_KTjet sub-leading jets pair invariant mass", 11, 50, 160); 

    TH1D *AKTjetMass1_bkg4 = new TH1D("AKTjetMass1_bkg4", "Anti_KTjet leading jets pair invariant mass", 11, 50, 160);
    TH1D *AKTjetMass2_bkg4 = new TH1D("AKTjetMass2_bkg4", "Anti_KTjet sub-leading jets pair invariant mass", 11, 50, 160); 
*/

    Pairing_w_JES(inputSigFile, JER, JERR02, AKTjetMass1_sig, AKTjetMass2_sig, tree_BDT_sig, output); 
    Pairing_w_JES(inputBkg1File, JER, JERR02, AKTjetMass1_bkg1, AKTjetMass2_bkg1, tree_BDT_bkg1, output); 
    Pairing_w_JES(inputBkg2File, JER, JERR02, AKTjetMass1_bkg2, AKTjetMass2_bkg2, tree_BDT_bkg2, output); 
    //Pairing_w_JES(inputBkg3File, JER, AKTjetMass1_bkg3, AKTjetMass2_bkg3, tree_BDT_bkg3, output); 
    //Pairing_w_JES(inputBkg4File, JER, AKTjetMass1_bkg4, AKTjetMass2_bkg4, tree_BDT_bkg4, output); 

    output -> cd();

    //JetPair1 -> Write();
    //JetPair2 -> Write();

    //tree_output -> Write();
    tree_BDT_sig -> Write();
    tree_BDT_bkg1 -> Write();
    tree_BDT_bkg2 -> Write();
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
