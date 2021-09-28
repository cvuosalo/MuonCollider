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

    //Calibration histogram
    TH1D *jetPTresponse = new TH1D("jetPTreponse", "jetPTresponse", 50, 0, 2);
    TH1D *MuonwjetPTresponse = new TH1D("MuonwjetPTreponse", "MuonwjetPTresponse", 50, 0, 2);

    TH1D *jetPTresponse00 = new TH1D("jetPTreponse00", "jetPTresponse00", 50, 0, 3.3);
    TH1D *jetPTresponse01 = new TH1D("jetPTreponse01", "jetPTresponse01", 50, 0, 3.3);
    TH1D *jetPTresponse02 = new TH1D("jetPTreponse02", "jetPTresponse02", 50, 0, 3.3);
    TH1D *jetPTresponse03 = new TH1D("jetPTreponse03", "jetPTresponse03", 50, 0, 3.3);
    TH1D *jetPTresponse04 = new TH1D("jetPTreponse04", "jetPTresponse04", 50, 0, 3.3);
    TH1D *jetPTresponse05 = new TH1D("jetPTreponse05", "jetPTresponse05", 50, 0, 3.3);
    TH1D *jetPTresponse06 = new TH1D("jetPTreponse06", "jetPTresponse06", 50, 0, 3.3);
    TH1D *jetPTresponse07 = new TH1D("jetPTreponse07", "jetPTresponse07", 50, 0, 3.3);
    TH1D *jetPTresponse08 = new TH1D("jetPTreponse08", "jetPTresponse08", 50, 0, 3.3);
    TH1D *jetPTresponse09 = new TH1D("jetPTreponse09", "jetPTresponse09", 50, 0, 3.3);

    TH1D *jetPTresponse10 = new TH1D("jetPTreponse10", "jetPTresponse10", 50, 0, 3.3);
    TH1D *jetPTresponse11 = new TH1D("jetPTreponse11", "jetPTresponse11", 50, 0, 3.3);
    TH1D *jetPTresponse12 = new TH1D("jetPTreponse12", "jetPTresponse12", 50, 0, 3.3);
    TH1D *jetPTresponse13 = new TH1D("jetPTreponse13", "jetPTresponse13", 50, 0, 3.3);
    TH1D *jetPTresponse14 = new TH1D("jetPTreponse14", "jetPTresponse14", 50, 0, 3.3);
    TH1D *jetPTresponse15 = new TH1D("jetPTreponse15", "jetPTresponse15", 50, 0, 3.3);
    TH1D *jetPTresponse16 = new TH1D("jetPTreponse16", "jetPTresponse16", 50, 0, 3.3);
    TH1D *jetPTresponse17 = new TH1D("jetPTreponse17", "jetPTresponse17", 50, 0, 3.3);
    TH1D *jetPTresponse18 = new TH1D("jetPTreponse18", "jetPTresponse18", 50, 0, 3.3);
    TH1D *jetPTresponse19 = new TH1D("jetPTreponse19", "jetPTresponse19", 50, 0, 3.3);

    TH1D *jetPTresponse20 = new TH1D("jetPTreponse20", "jetPTresponse20", 50, 0, 3.3);
    TH1D *jetPTresponse21 = new TH1D("jetPTreponse21", "jetPTresponse21", 50, 0, 3.3);
    TH1D *jetPTresponse22 = new TH1D("jetPTreponse22", "jetPTresponse22", 50, 0, 3.3);
    TH1D *jetPTresponse23 = new TH1D("jetPTreponse23", "jetPTresponse23", 50, 0, 3.3);
    TH1D *jetPTresponse24 = new TH1D("jetPTreponse24", "jetPTresponse24", 50, 0, 3.3);
    TH1D *jetPTresponse25 = new TH1D("jetPTreponse25", "jetPTresponse25", 50, 0, 3.3);
    TH1D *jetPTresponse26 = new TH1D("jetPTreponse26", "jetPTresponse26", 50, 0, 3.3);
    TH1D *jetPTresponse27 = new TH1D("jetPTreponse27", "jetPTresponse27", 50, 0, 3.3);
    TH1D *jetPTresponse28 = new TH1D("jetPTreponse28", "jetPTresponse28", 50, 0, 3.3);
    TH1D *jetPTresponse29 = new TH1D("jetPTreponse29", "jetPTresponse29", 50, 0, 3.3);

    TH1D *jetPTresponse30 = new TH1D("jetPTreponse30", "jetPTresponse30", 50, 0, 3.3);
    TH1D *jetPTresponse31 = new TH1D("jetPTreponse31", "jetPTresponse31", 50, 0, 3.3);
    TH1D *jetPTresponse32 = new TH1D("jetPTreponse32", "jetPTresponse32", 50, 0, 3.3);
    TH1D *jetPTresponse33 = new TH1D("jetPTreponse33", "jetPTresponse33", 50, 0, 3.3);
    TH1D *jetPTresponse34 = new TH1D("jetPTreponse34", "jetPTresponse34", 50, 0, 3.3);
    TH1D *jetPTresponse35 = new TH1D("jetPTreponse35", "jetPTresponse35", 50, 0, 3.3);
    TH1D *jetPTresponse36 = new TH1D("jetPTreponse36", "jetPTresponse36", 50, 0, 3.3);
    TH1D *jetPTresponse37 = new TH1D("jetPTreponse37", "jetPTresponse37", 50, 0, 3.3);
    TH1D *jetPTresponse38 = new TH1D("jetPTreponse38", "jetPTresponse38", 50, 0, 3.3);
    TH1D *jetPTresponse39 = new TH1D("jetPTreponse39", "jetPTresponse39", 50, 0, 3.3);

    TH1D *jetPTresponse40 = new TH1D("jetPTreponse40", "jetPTresponse40", 50, 0, 3.3);
    TH1D *jetPTresponse41 = new TH1D("jetPTreponse41", "jetPTresponse41", 50, 0, 3.3);
    TH1D *jetPTresponse42 = new TH1D("jetPTreponse42", "jetPTresponse42", 50, 0, 3.3);
    TH1D *jetPTresponse43 = new TH1D("jetPTreponse43", "jetPTresponse43", 50, 0, 3.3);
    TH1D *jetPTresponse44 = new TH1D("jetPTreponse44", "jetPTresponse44", 50, 0, 3.3);
    TH1D *jetPTresponse45 = new TH1D("jetPTreponse45", "jetPTresponse45", 50, 0, 3.3);
    TH1D *jetPTresponse46 = new TH1D("jetPTreponse46", "jetPTresponse46", 50, 0, 3.3);
    TH1D *jetPTresponse47 = new TH1D("jetPTreponse47", "jetPTresponse47", 50, 0, 3.3);
    TH1D *jetPTresponse48 = new TH1D("jetPTreponse48", "jetPTresponse48", 50, 0, 3.3);
    TH1D *jetPTresponse49 = new TH1D("jetPTreponse49", "jetPTresponse49", 50, 0, 3.3);

    TH1D *jetPTresponse50 = new TH1D("jetPTreponse50", "jetPTresponse50", 50, 0, 3.3);
    TH1D *jetPTresponse51 = new TH1D("jetPTreponse51", "jetPTresponse51", 50, 0, 3.3);
    TH1D *jetPTresponse52 = new TH1D("jetPTreponse52", "jetPTresponse52", 50, 0, 3.3);
    TH1D *jetPTresponse53 = new TH1D("jetPTreponse53", "jetPTresponse53", 50, 0, 3.3);
    TH1D *jetPTresponse54 = new TH1D("jetPTreponse54", "jetPTresponse54", 50, 0, 3.3);
    TH1D *jetPTresponse55 = new TH1D("jetPTreponse55", "jetPTresponse55", 50, 0, 3.3);
    TH1D *jetPTresponse56 = new TH1D("jetPTreponse56", "jetPTresponse56", 50, 0, 3.3);
    TH1D *jetPTresponse57 = new TH1D("jetPTreponse57", "jetPTresponse57", 50, 0, 3.3);
    TH1D *jetPTresponse58 = new TH1D("jetPTreponse58", "jetPTresponse58", 50, 0, 3.3);
    TH1D *jetPTresponse59 = new TH1D("jetPTreponse59", "jetPTresponse59", 50, 0, 3.3);

    TH1D *jetPTresponse60 = new TH1D("jetPTreponse60", "jetPTresponse60", 50, 0, 3.3);
    TH1D *jetPTresponse61 = new TH1D("jetPTreponse61", "jetPTresponse61", 50, 0, 3.3);
    TH1D *jetPTresponse62 = new TH1D("jetPTreponse62", "jetPTresponse62", 50, 0, 3.3);
    TH1D *jetPTresponse63 = new TH1D("jetPTreponse63", "jetPTresponse63", 50, 0, 3.3);
    TH1D *jetPTresponse64 = new TH1D("jetPTreponse64", "jetPTresponse64", 50, 0, 3.3);
    TH1D *jetPTresponse65 = new TH1D("jetPTreponse65", "jetPTresponse65", 50, 0, 3.3);
    TH1D *jetPTresponse66 = new TH1D("jetPTreponse66", "jetPTresponse66", 50, 0, 3.3);
    TH1D *jetPTresponse67 = new TH1D("jetPTreponse67", "jetPTresponse67", 50, 0, 3.3);
    TH1D *jetPTresponse68 = new TH1D("jetPTreponse68", "jetPTresponse68", 50, 0, 3.3);
    TH1D *jetPTresponse69 = new TH1D("jetPTreponse69", "jetPTresponse69", 50, 0, 3.3);

    TH1D *jetPTresponse70 = new TH1D("jetPTreponse70", "jetPTresponse70", 50, 0, 3.3);
    TH1D *jetPTresponse71 = new TH1D("jetPTreponse71", "jetPTresponse71", 50, 0, 3.3);
    TH1D *jetPTresponse72 = new TH1D("jetPTreponse72", "jetPTresponse72", 50, 0, 3.3);
    TH1D *jetPTresponse73 = new TH1D("jetPTreponse73", "jetPTresponse73", 50, 0, 3.3);
    TH1D *jetPTresponse74 = new TH1D("jetPTreponse74", "jetPTresponse74", 50, 0, 3.3);
    TH1D *jetPTresponse75 = new TH1D("jetPTreponse75", "jetPTresponse75", 50, 0, 3.3);
    TH1D *jetPTresponse76 = new TH1D("jetPTreponse76", "jetPTresponse76", 50, 0, 3.3);
    TH1D *jetPTresponse77 = new TH1D("jetPTreponse77", "jetPTresponse77", 50, 0, 3.3);
    TH1D *jetPTresponse78 = new TH1D("jetPTreponse78", "jetPTresponse78", 50, 0, 3.3);
    TH1D *jetPTresponse79 = new TH1D("jetPTreponse79", "jetPTresponse79", 50, 0, 3.3);

    TH1D *jetPTresponse80 = new TH1D("jetPTreponse80", "jetPTresponse80", 50, 0, 3.3);
    TH1D *jetPTresponse81 = new TH1D("jetPTreponse81", "jetPTresponse81", 50, 0, 3.3);
    TH1D *jetPTresponse82 = new TH1D("jetPTreponse82", "jetPTresponse82", 50, 0, 3.3);
    TH1D *jetPTresponse83 = new TH1D("jetPTreponse83", "jetPTresponse83", 50, 0, 3.3);
    TH1D *jetPTresponse84 = new TH1D("jetPTreponse84", "jetPTresponse84", 50, 0, 3.3);
    TH1D *jetPTresponse85 = new TH1D("jetPTreponse85", "jetPTresponse85", 50, 0, 3.3);
    TH1D *jetPTresponse86 = new TH1D("jetPTreponse86", "jetPTresponse86", 50, 0, 3.3);
    TH1D *jetPTresponse87 = new TH1D("jetPTreponse87", "jetPTresponse87", 50, 0, 3.3);
    TH1D *jetPTresponse88 = new TH1D("jetPTreponse88", "jetPTresponse88", 50, 0, 3.3);
    TH1D *jetPTresponse89 = new TH1D("jetPTreponse89", "jetPTresponse89", 50, 0, 3.3);

    TH1D *jetPTresponse90 = new TH1D("jetPTreponse90", "jetPTresponse90", 50, 0, 3.3);
    TH1D *jetPTresponse91 = new TH1D("jetPTreponse91", "jetPTresponse91", 50, 0, 3.3);
    TH1D *jetPTresponse92 = new TH1D("jetPTreponse92", "jetPTresponse92", 50, 0, 3.3);
    TH1D *jetPTresponse93 = new TH1D("jetPTreponse93", "jetPTresponse93", 50, 0, 3.3);
    TH1D *jetPTresponse94 = new TH1D("jetPTreponse94", "jetPTresponse94", 50, 0, 3.3);
    TH1D *jetPTresponse95 = new TH1D("jetPTreponse95", "jetPTresponse95", 50, 0, 3.3);
    TH1D *jetPTresponse96 = new TH1D("jetPTreponse96", "jetPTresponse96", 50, 0, 3.3);
    TH1D *jetPTresponse97 = new TH1D("jetPTreponse97", "jetPTresponse97", 50, 0, 3.3);
    TH1D *jetPTresponse98 = new TH1D("jetPTreponse98", "jetPTresponse98", 50, 0, 3.3);
    TH1D *jetPTresponse99 = new TH1D("jetPTreponse99", "jetPTresponse99", 50, 0, 3.3);

    TH1D *MuonwjetPTresponse00 = new TH1D("MuonwjetPTreponse00", "MuonwjetPTresponse00", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse01 = new TH1D("MuonwjetPTreponse01", "MuonwjetPTresponse01", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse02 = new TH1D("MuonwjetPTreponse02", "MuonwjetPTresponse02", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse03 = new TH1D("MuonwjetPTreponse03", "MuonwjetPTresponse03", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse04 = new TH1D("MuonwjetPTreponse04", "MuonwjetPTresponse04", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse05 = new TH1D("MuonwjetPTreponse05", "MuonwjetPTresponse05", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse06 = new TH1D("MuonwjetPTreponse06", "MuonwjetPTresponse06", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse07 = new TH1D("MuonwjetPTreponse07", "MuonwjetPTresponse07", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse08 = new TH1D("MuonwjetPTreponse08", "MuonwjetPTresponse08", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse09 = new TH1D("MuonwjetPTreponse09", "MuonwjetPTresponse09", 50, 0, 3.3);

    TH1D *MuonwjetPTresponse10 = new TH1D("MuonwjetPTreponse10", "MuonwjetPTresponse10", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse11 = new TH1D("MuonwjetPTreponse11", "MuonwjetPTresponse11", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse12 = new TH1D("MuonwjetPTreponse12", "MuonwjetPTresponse12", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse13 = new TH1D("MuonwjetPTreponse13", "MuonwjetPTresponse13", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse14 = new TH1D("MuonwjetPTreponse14", "MuonwjetPTresponse14", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse15 = new TH1D("MuonwjetPTreponse15", "MuonwjetPTresponse15", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse16 = new TH1D("MuonwjetPTreponse16", "MuonwjetPTresponse16", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse17 = new TH1D("MuonwjetPTreponse17", "MuonwjetPTresponse17", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse18 = new TH1D("MuonwjetPTreponse18", "MuonwjetPTresponse18", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse19 = new TH1D("MuonwjetPTreponse19", "MuonwjetPTresponse19", 50, 0, 3.3);

    TH1D *MuonwjetPTresponse20 = new TH1D("MuonwjetPTreponse20", "MuonwjetPTresponse20", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse21 = new TH1D("MuonwjetPTreponse21", "MuonwjetPTresponse21", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse22 = new TH1D("MuonwjetPTreponse22", "MuonwjetPTresponse22", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse23 = new TH1D("MuonwjetPTreponse23", "MuonwjetPTresponse23", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse24 = new TH1D("MuonwjetPTreponse24", "MuonwjetPTresponse24", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse25 = new TH1D("MuonwjetPTreponse25", "MuonwjetPTresponse25", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse26 = new TH1D("MuonwjetPTreponse26", "MuonwjetPTresponse26", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse27 = new TH1D("MuonwjetPTreponse27", "MuonwjetPTresponse27", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse28 = new TH1D("MuonwjetPTreponse28", "MuonwjetPTresponse28", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse29 = new TH1D("MuonwjetPTreponse29", "MuonwjetPTresponse29", 50, 0, 3.3);

    TH1D *MuonwjetPTresponse30 = new TH1D("MuonwjetPTreponse30", "MuonwjetPTresponse30", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse31 = new TH1D("MuonwjetPTreponse31", "MuonwjetPTresponse31", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse32 = new TH1D("MuonwjetPTreponse32", "MuonwjetPTresponse32", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse33 = new TH1D("MuonwjetPTreponse33", "MuonwjetPTresponse33", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse34 = new TH1D("MuonwjetPTreponse34", "MuonwjetPTresponse34", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse35 = new TH1D("MuonwjetPTreponse35", "MuonwjetPTresponse35", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse36 = new TH1D("MuonwjetPTreponse36", "MuonwjetPTresponse36", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse37 = new TH1D("MuonwjetPTreponse37", "MuonwjetPTresponse37", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse38 = new TH1D("MuonwjetPTreponse38", "MuonwjetPTresponse38", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse39 = new TH1D("MuonwjetPTreponse39", "MuonwjetPTresponse39", 50, 0, 3.3);

    TH1D *MuonwjetPTresponse40 = new TH1D("MuonwjetPTreponse40", "MuonwjetPTresponse40", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse41 = new TH1D("MuonwjetPTreponse41", "MuonwjetPTresponse41", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse42 = new TH1D("MuonwjetPTreponse42", "MuonwjetPTresponse42", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse43 = new TH1D("MuonwjetPTreponse43", "MuonwjetPTresponse43", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse44 = new TH1D("MuonwjetPTreponse44", "MuonwjetPTresponse44", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse45 = new TH1D("MuonwjetPTreponse45", "MuonwjetPTresponse45", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse46 = new TH1D("MuonwjetPTreponse46", "MuonwjetPTresponse46", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse47 = new TH1D("MuonwjetPTreponse47", "MuonwjetPTresponse47", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse48 = new TH1D("MuonwjetPTreponse48", "MuonwjetPTresponse48", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse49 = new TH1D("MuonwjetPTreponse49", "MuonwjetPTresponse49", 50, 0, 3.3);

    TH1D *MuonwjetPTresponse50 = new TH1D("MuonwjetPTreponse50", "MuonwjetPTresponse50", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse51 = new TH1D("MuonwjetPTreponse51", "MuonwjetPTresponse51", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse52 = new TH1D("MuonwjetPTreponse52", "MuonwjetPTresponse52", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse53 = new TH1D("MuonwjetPTreponse53", "MuonwjetPTresponse53", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse54 = new TH1D("MuonwjetPTreponse54", "MuonwjetPTresponse54", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse55 = new TH1D("MuonwjetPTreponse55", "MuonwjetPTresponse55", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse56 = new TH1D("MuonwjetPTreponse56", "MuonwjetPTresponse56", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse57 = new TH1D("MuonwjetPTreponse57", "MuonwjetPTresponse57", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse58 = new TH1D("MuonwjetPTreponse58", "MuonwjetPTresponse58", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse59 = new TH1D("MuonwjetPTreponse59", "MuonwjetPTresponse59", 50, 0, 3.3);

    TH1D *MuonwjetPTresponse60 = new TH1D("MuonwjetPTreponse60", "MuonwjetPTresponse60", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse61 = new TH1D("MuonwjetPTreponse61", "MuonwjetPTresponse61", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse62 = new TH1D("MuonwjetPTreponse62", "MuonwjetPTresponse62", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse63 = new TH1D("MuonwjetPTreponse63", "MuonwjetPTresponse63", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse64 = new TH1D("MuonwjetPTreponse64", "MuonwjetPTresponse64", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse65 = new TH1D("MuonwjetPTreponse65", "MuonwjetPTresponse65", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse66 = new TH1D("MuonwjetPTreponse66", "MuonwjetPTresponse66", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse67 = new TH1D("MuonwjetPTreponse67", "MuonwjetPTresponse67", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse68 = new TH1D("MuonwjetPTreponse68", "MuonwjetPTresponse68", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse69 = new TH1D("MuonwjetPTreponse69", "MuonwjetPTresponse69", 50, 0, 3.3);

    TH1D *MuonwjetPTresponse70 = new TH1D("MuonwjetPTreponse70", "MuonwjetPTresponse70", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse71 = new TH1D("MuonwjetPTreponse71", "MuonwjetPTresponse71", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse72 = new TH1D("MuonwjetPTreponse72", "MuonwjetPTresponse72", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse73 = new TH1D("MuonwjetPTreponse73", "MuonwjetPTresponse73", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse74 = new TH1D("MuonwjetPTreponse74", "MuonwjetPTresponse74", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse75 = new TH1D("MuonwjetPTreponse75", "MuonwjetPTresponse75", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse76 = new TH1D("MuonwjetPTreponse76", "MuonwjetPTresponse76", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse77 = new TH1D("MuonwjetPTreponse77", "MuonwjetPTresponse77", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse78 = new TH1D("MuonwjetPTreponse78", "MuonwjetPTresponse78", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse79 = new TH1D("MuonwjetPTreponse79", "MuonwjetPTresponse79", 50, 0, 3.3);

    TH1D *MuonwjetPTresponse80 = new TH1D("MuonwjetPTreponse80", "MuonwjetPTresponse80", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse81 = new TH1D("MuonwjetPTreponse81", "MuonwjetPTresponse81", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse82 = new TH1D("MuonwjetPTreponse82", "MuonwjetPTresponse82", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse83 = new TH1D("MuonwjetPTreponse83", "MuonwjetPTresponse83", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse84 = new TH1D("MuonwjetPTreponse84", "MuonwjetPTresponse84", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse85 = new TH1D("MuonwjetPTreponse85", "MuonwjetPTresponse85", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse86 = new TH1D("MuonwjetPTreponse86", "MuonwjetPTresponse86", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse87 = new TH1D("MuonwjetPTreponse87", "MuonwjetPTresponse87", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse88 = new TH1D("MuonwjetPTreponse88", "MuonwjetPTresponse88", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse89 = new TH1D("MuonwjetPTreponse89", "MuonwjetPTresponse89", 50, 0, 3.3);

    TH1D *MuonwjetPTresponse90 = new TH1D("MuonwjetPTreponse90", "MuonwjetPTresponse90", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse91 = new TH1D("MuonwjetPTreponse91", "MuonwjetPTresponse91", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse92 = new TH1D("MuonwjetPTreponse92", "MuonwjetPTresponse92", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse93 = new TH1D("MuonwjetPTreponse93", "MuonwjetPTresponse93", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse94 = new TH1D("MuonwjetPTreponse94", "MuonwjetPTresponse94", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse95 = new TH1D("MuonwjetPTreponse95", "MuonwjetPTresponse95", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse96 = new TH1D("MuonwjetPTreponse96", "MuonwjetPTresponse96", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse97 = new TH1D("MuonwjetPTreponse97", "MuonwjetPTresponse97", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse98 = new TH1D("MuonwjetPTreponse98", "MuonwjetPTresponse98", 50, 0, 3.3);
    TH1D *MuonwjetPTresponse99 = new TH1D("MuonwjetPTreponse99", "MuonwjetPTresponse99", 50, 0, 3.3);

    TH1D *MuonwojetPTresponse00 = new TH1D("MuonwojetPTreponse00", "MuonwojetPTresponse00", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse01 = new TH1D("MuonwojetPTreponse01", "MuonwojetPTresponse01", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse02 = new TH1D("MuonwojetPTreponse02", "MuonwojetPTresponse02", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse03 = new TH1D("MuonwojetPTreponse03", "MuonwojetPTresponse03", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse04 = new TH1D("MuonwojetPTreponse04", "MuonwojetPTresponse04", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse05 = new TH1D("MuonwojetPTreponse05", "MuonwojetPTresponse05", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse06 = new TH1D("MuonwojetPTreponse06", "MuonwojetPTresponse06", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse07 = new TH1D("MuonwojetPTreponse07", "MuonwojetPTresponse07", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse08 = new TH1D("MuonwojetPTreponse08", "MuonwojetPTresponse08", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse09 = new TH1D("MuonwojetPTreponse09", "MuonwojetPTresponse09", 50, 0, 3.3);

    TH1D *MuonwojetPTresponse10 = new TH1D("MuonwojetPTreponse10", "MuonwojetPTresponse10", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse11 = new TH1D("MuonwojetPTreponse11", "MuonwojetPTresponse11", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse12 = new TH1D("MuonwojetPTreponse12", "MuonwojetPTresponse12", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse13 = new TH1D("MuonwojetPTreponse13", "MuonwojetPTresponse13", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse14 = new TH1D("MuonwojetPTreponse14", "MuonwojetPTresponse14", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse15 = new TH1D("MuonwojetPTreponse15", "MuonwojetPTresponse15", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse16 = new TH1D("MuonwojetPTreponse16", "MuonwojetPTresponse16", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse17 = new TH1D("MuonwojetPTreponse17", "MuonwojetPTresponse17", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse18 = new TH1D("MuonwojetPTreponse18", "MuonwojetPTresponse18", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse19 = new TH1D("MuonwojetPTreponse19", "MuonwojetPTresponse19", 50, 0, 3.3);

    TH1D *MuonwojetPTresponse20 = new TH1D("MuonwojetPTreponse20", "MuonwojetPTresponse20", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse21 = new TH1D("MuonwojetPTreponse21", "MuonwojetPTresponse21", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse22 = new TH1D("MuonwojetPTreponse22", "MuonwojetPTresponse22", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse23 = new TH1D("MuonwojetPTreponse23", "MuonwojetPTresponse23", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse24 = new TH1D("MuonwojetPTreponse24", "MuonwojetPTresponse24", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse25 = new TH1D("MuonwojetPTreponse25", "MuonwojetPTresponse25", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse26 = new TH1D("MuonwojetPTreponse26", "MuonwojetPTresponse26", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse27 = new TH1D("MuonwojetPTreponse27", "MuonwojetPTresponse27", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse28 = new TH1D("MuonwojetPTreponse28", "MuonwojetPTresponse28", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse29 = new TH1D("MuonwojetPTreponse29", "MuonwojetPTresponse29", 50, 0, 3.3);

    TH1D *MuonwojetPTresponse30 = new TH1D("MuonwojetPTreponse30", "MuonwojetPTresponse30", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse31 = new TH1D("MuonwojetPTreponse31", "MuonwojetPTresponse31", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse32 = new TH1D("MuonwojetPTreponse32", "MuonwojetPTresponse32", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse33 = new TH1D("MuonwojetPTreponse33", "MuonwojetPTresponse33", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse34 = new TH1D("MuonwojetPTreponse34", "MuonwojetPTresponse34", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse35 = new TH1D("MuonwojetPTreponse35", "MuonwojetPTresponse35", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse36 = new TH1D("MuonwojetPTreponse36", "MuonwojetPTresponse36", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse37 = new TH1D("MuonwojetPTreponse37", "MuonwojetPTresponse37", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse38 = new TH1D("MuonwojetPTreponse38", "MuonwojetPTresponse38", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse39 = new TH1D("MuonwojetPTreponse39", "MuonwojetPTresponse39", 50, 0, 3.3);

    TH1D *MuonwojetPTresponse40 = new TH1D("MuonwojetPTreponse40", "MuonwojetPTresponse40", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse41 = new TH1D("MuonwojetPTreponse41", "MuonwojetPTresponse41", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse42 = new TH1D("MuonwojetPTreponse42", "MuonwojetPTresponse42", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse43 = new TH1D("MuonwojetPTreponse43", "MuonwojetPTresponse43", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse44 = new TH1D("MuonwojetPTreponse44", "MuonwojetPTresponse44", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse45 = new TH1D("MuonwojetPTreponse45", "MuonwojetPTresponse45", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse46 = new TH1D("MuonwojetPTreponse46", "MuonwojetPTresponse46", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse47 = new TH1D("MuonwojetPTreponse47", "MuonwojetPTresponse47", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse48 = new TH1D("MuonwojetPTreponse48", "MuonwojetPTresponse48", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse49 = new TH1D("MuonwojetPTreponse49", "MuonwojetPTresponse49", 50, 0, 3.3);

    TH1D *MuonwojetPTresponse50 = new TH1D("MuonwojetPTreponse50", "MuonwojetPTresponse50", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse51 = new TH1D("MuonwojetPTreponse51", "MuonwojetPTresponse51", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse52 = new TH1D("MuonwojetPTreponse52", "MuonwojetPTresponse52", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse53 = new TH1D("MuonwojetPTreponse53", "MuonwojetPTresponse53", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse54 = new TH1D("MuonwojetPTreponse54", "MuonwojetPTresponse54", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse55 = new TH1D("MuonwojetPTreponse55", "MuonwojetPTresponse55", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse56 = new TH1D("MuonwojetPTreponse56", "MuonwojetPTresponse56", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse57 = new TH1D("MuonwojetPTreponse57", "MuonwojetPTresponse57", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse58 = new TH1D("MuonwojetPTreponse58", "MuonwojetPTresponse58", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse59 = new TH1D("MuonwojetPTreponse59", "MuonwojetPTresponse59", 50, 0, 3.3);

    TH1D *MuonwojetPTresponse60 = new TH1D("MuonwojetPTreponse60", "MuonwojetPTresponse60", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse61 = new TH1D("MuonwojetPTreponse61", "MuonwojetPTresponse61", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse62 = new TH1D("MuonwojetPTreponse62", "MuonwojetPTresponse62", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse63 = new TH1D("MuonwojetPTreponse63", "MuonwojetPTresponse63", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse64 = new TH1D("MuonwojetPTreponse64", "MuonwojetPTresponse64", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse65 = new TH1D("MuonwojetPTreponse65", "MuonwojetPTresponse65", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse66 = new TH1D("MuonwojetPTreponse66", "MuonwojetPTresponse66", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse67 = new TH1D("MuonwojetPTreponse67", "MuonwojetPTresponse67", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse68 = new TH1D("MuonwojetPTreponse68", "MuonwojetPTresponse68", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse69 = new TH1D("MuonwojetPTreponse69", "MuonwojetPTresponse69", 50, 0, 3.3);

    TH1D *MuonwojetPTresponse70 = new TH1D("MuonwojetPTreponse70", "MuonwojetPTresponse70", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse71 = new TH1D("MuonwojetPTreponse71", "MuonwojetPTresponse71", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse72 = new TH1D("MuonwojetPTreponse72", "MuonwojetPTresponse72", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse73 = new TH1D("MuonwojetPTreponse73", "MuonwojetPTresponse73", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse74 = new TH1D("MuonwojetPTreponse74", "MuonwojetPTresponse74", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse75 = new TH1D("MuonwojetPTreponse75", "MuonwojetPTresponse75", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse76 = new TH1D("MuonwojetPTreponse76", "MuonwojetPTresponse76", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse77 = new TH1D("MuonwojetPTreponse77", "MuonwojetPTresponse77", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse78 = new TH1D("MuonwojetPTreponse78", "MuonwojetPTresponse78", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse79 = new TH1D("MuonwojetPTreponse79", "MuonwojetPTresponse79", 50, 0, 3.3);

    TH1D *MuonwojetPTresponse80 = new TH1D("MuonwojetPTreponse80", "MuonwojetPTresponse80", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse81 = new TH1D("MuonwojetPTreponse81", "MuonwojetPTresponse81", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse82 = new TH1D("MuonwojetPTreponse82", "MuonwojetPTresponse82", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse83 = new TH1D("MuonwojetPTreponse83", "MuonwojetPTresponse83", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse84 = new TH1D("MuonwojetPTreponse84", "MuonwojetPTresponse84", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse85 = new TH1D("MuonwojetPTreponse85", "MuonwojetPTresponse85", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse86 = new TH1D("MuonwojetPTreponse86", "MuonwojetPTresponse86", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse87 = new TH1D("MuonwojetPTreponse87", "MuonwojetPTresponse87", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse88 = new TH1D("MuonwojetPTreponse88", "MuonwojetPTresponse88", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse89 = new TH1D("MuonwojetPTreponse89", "MuonwojetPTresponse89", 50, 0, 3.3);

    TH1D *MuonwojetPTresponse90 = new TH1D("MuonwojetPTreponse90", "MuonwojetPTresponse90", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse91 = new TH1D("MuonwojetPTreponse91", "MuonwojetPTresponse91", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse92 = new TH1D("MuonwojetPTreponse92", "MuonwojetPTresponse92", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse93 = new TH1D("MuonwojetPTreponse93", "MuonwojetPTresponse93", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse94 = new TH1D("MuonwojetPTreponse94", "MuonwojetPTresponse94", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse95 = new TH1D("MuonwojetPTreponse95", "MuonwojetPTresponse95", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse96 = new TH1D("MuonwojetPTreponse96", "MuonwojetPTresponse96", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse97 = new TH1D("MuonwojetPTreponse97", "MuonwojetPTresponse97", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse98 = new TH1D("MuonwojetPTreponse98", "MuonwojetPTresponse98", 50, 0, 3.3);
    TH1D *MuonwojetPTresponse99 = new TH1D("MuonwojetPTreponse99", "MuonwojetPTresponse99", 50, 0, 3.3);

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
    TH1D *jetPT2darray[10][10][3];
    bool MuonTagging = false;

    //Calibration 2darray
    jetPT2darray[0][0][0] = jetPTresponse00;
    jetPT2darray[0][1][0] = jetPTresponse01;
    jetPT2darray[0][2][0] = jetPTresponse02;
    jetPT2darray[0][3][0] = jetPTresponse03;
    jetPT2darray[0][4][0] = jetPTresponse04;
    jetPT2darray[0][5][0] = jetPTresponse05;
    jetPT2darray[0][6][0] = jetPTresponse06;
    jetPT2darray[0][7][0] = jetPTresponse07;
    jetPT2darray[0][8][0] = jetPTresponse08;
    jetPT2darray[0][9][0] = jetPTresponse09;

    jetPT2darray[1][0][0] = jetPTresponse10;
    jetPT2darray[1][1][0] = jetPTresponse11;
    jetPT2darray[1][2][0] = jetPTresponse12;
    jetPT2darray[1][3][0] = jetPTresponse13;
    jetPT2darray[1][4][0] = jetPTresponse14;
    jetPT2darray[1][5][0] = jetPTresponse15;
    jetPT2darray[1][6][0] = jetPTresponse16;
    jetPT2darray[1][7][0] = jetPTresponse17;
    jetPT2darray[1][8][0] = jetPTresponse18;
    jetPT2darray[1][9][0] = jetPTresponse19;

    jetPT2darray[2][0][0] = jetPTresponse20;
    jetPT2darray[2][1][0] = jetPTresponse21;
    jetPT2darray[2][2][0] = jetPTresponse22;
    jetPT2darray[2][3][0] = jetPTresponse23;
    jetPT2darray[2][4][0] = jetPTresponse24;
    jetPT2darray[2][5][0] = jetPTresponse25;
    jetPT2darray[2][6][0] = jetPTresponse26;
    jetPT2darray[2][7][0] = jetPTresponse27;
    jetPT2darray[2][8][0] = jetPTresponse28;
    jetPT2darray[2][9][0] = jetPTresponse29;

    jetPT2darray[3][0][0] = jetPTresponse30;
    jetPT2darray[3][1][0] = jetPTresponse31;
    jetPT2darray[3][2][0] = jetPTresponse32;
    jetPT2darray[3][3][0] = jetPTresponse33;
    jetPT2darray[3][4][0] = jetPTresponse34;
    jetPT2darray[3][5][0] = jetPTresponse35;
    jetPT2darray[3][6][0] = jetPTresponse36;
    jetPT2darray[3][7][0] = jetPTresponse37;
    jetPT2darray[3][8][0] = jetPTresponse38;
    jetPT2darray[3][9][0] = jetPTresponse39;

    jetPT2darray[4][0][0] = jetPTresponse40;
    jetPT2darray[4][1][0] = jetPTresponse41;
    jetPT2darray[4][2][0] = jetPTresponse42;
    jetPT2darray[4][3][0] = jetPTresponse43;
    jetPT2darray[4][4][0] = jetPTresponse44;
    jetPT2darray[4][5][0] = jetPTresponse45;
    jetPT2darray[4][6][0] = jetPTresponse46;
    jetPT2darray[4][7][0] = jetPTresponse47;
    jetPT2darray[4][8][0] = jetPTresponse48;
    jetPT2darray[4][9][0] = jetPTresponse49;

    jetPT2darray[5][0][0] = jetPTresponse50;
    jetPT2darray[5][1][0] = jetPTresponse51;
    jetPT2darray[5][2][0] = jetPTresponse52;
    jetPT2darray[5][3][0] = jetPTresponse53;
    jetPT2darray[5][4][0] = jetPTresponse54;
    jetPT2darray[5][5][0] = jetPTresponse55;
    jetPT2darray[5][6][0] = jetPTresponse56;
    jetPT2darray[5][7][0] = jetPTresponse57;
    jetPT2darray[5][8][0] = jetPTresponse58;
    jetPT2darray[5][9][0] = jetPTresponse59;

    jetPT2darray[6][0][0] = jetPTresponse60;
    jetPT2darray[6][1][0] = jetPTresponse61;
    jetPT2darray[6][2][0] = jetPTresponse62;
    jetPT2darray[6][3][0] = jetPTresponse63;
    jetPT2darray[6][4][0] = jetPTresponse64;
    jetPT2darray[6][5][0] = jetPTresponse65;
    jetPT2darray[6][6][0] = jetPTresponse66;
    jetPT2darray[6][7][0] = jetPTresponse67;
    jetPT2darray[6][8][0] = jetPTresponse68;
    jetPT2darray[6][9][0] = jetPTresponse69;

    jetPT2darray[7][0][0] = jetPTresponse70;
    jetPT2darray[7][1][0] = jetPTresponse71;
    jetPT2darray[7][2][0] = jetPTresponse72;
    jetPT2darray[7][3][0] = jetPTresponse73;
    jetPT2darray[7][4][0] = jetPTresponse74;
    jetPT2darray[7][5][0] = jetPTresponse75;
    jetPT2darray[7][6][0] = jetPTresponse76;
    jetPT2darray[7][7][0] = jetPTresponse77;
    jetPT2darray[7][8][0] = jetPTresponse78;
    jetPT2darray[7][9][0] = jetPTresponse79;

    jetPT2darray[8][0][0] = jetPTresponse80;
    jetPT2darray[8][1][0] = jetPTresponse81;
    jetPT2darray[8][2][0] = jetPTresponse82;
    jetPT2darray[8][3][0] = jetPTresponse83;
    jetPT2darray[8][4][0] = jetPTresponse84;
    jetPT2darray[8][5][0] = jetPTresponse85;
    jetPT2darray[8][6][0] = jetPTresponse86;
    jetPT2darray[8][7][0] = jetPTresponse87;
    jetPT2darray[8][8][0] = jetPTresponse88;
    jetPT2darray[8][9][0] = jetPTresponse89;

    jetPT2darray[9][0][0] = jetPTresponse90;
    jetPT2darray[9][1][0] = jetPTresponse91;
    jetPT2darray[9][2][0] = jetPTresponse92;
    jetPT2darray[9][3][0] = jetPTresponse93;
    jetPT2darray[9][4][0] = jetPTresponse94;
    jetPT2darray[9][5][0] = jetPTresponse95;
    jetPT2darray[9][6][0] = jetPTresponse96;
    jetPT2darray[9][7][0] = jetPTresponse97;
    jetPT2darray[9][8][0] = jetPTresponse98;
    jetPT2darray[9][9][0] = jetPTresponse99;

    jetPT2darray[0][0][1] = MuonwjetPTresponse00;
    jetPT2darray[0][1][1] = MuonwjetPTresponse01;
    jetPT2darray[0][2][1] = MuonwjetPTresponse02;
    jetPT2darray[0][3][1] = MuonwjetPTresponse03;
    jetPT2darray[0][4][1] = MuonwjetPTresponse04;
    jetPT2darray[0][5][1] = MuonwjetPTresponse05;
    jetPT2darray[0][6][1] = MuonwjetPTresponse06;
    jetPT2darray[0][7][1] = MuonwjetPTresponse07;
    jetPT2darray[0][8][1] = MuonwjetPTresponse08;
    jetPT2darray[0][9][1] = MuonwjetPTresponse09;

    jetPT2darray[1][0][1] = MuonwjetPTresponse10;
    jetPT2darray[1][1][1] = MuonwjetPTresponse11;
    jetPT2darray[1][2][1] = MuonwjetPTresponse12;
    jetPT2darray[1][3][1] = MuonwjetPTresponse13;
    jetPT2darray[1][4][1] = MuonwjetPTresponse14;
    jetPT2darray[1][5][1] = MuonwjetPTresponse15;
    jetPT2darray[1][6][1] = MuonwjetPTresponse16;
    jetPT2darray[1][7][1] = MuonwjetPTresponse17;
    jetPT2darray[1][8][1] = MuonwjetPTresponse18;
    jetPT2darray[1][9][1] = MuonwjetPTresponse19;

    jetPT2darray[2][0][1] = MuonwjetPTresponse20;
    jetPT2darray[2][1][1] = MuonwjetPTresponse21;
    jetPT2darray[2][2][1] = MuonwjetPTresponse22;
    jetPT2darray[2][3][1] = MuonwjetPTresponse23;
    jetPT2darray[2][4][1] = MuonwjetPTresponse24;
    jetPT2darray[2][5][1] = MuonwjetPTresponse25;
    jetPT2darray[2][6][1] = MuonwjetPTresponse26;
    jetPT2darray[2][7][1] = MuonwjetPTresponse27;
    jetPT2darray[2][8][1] = MuonwjetPTresponse28;
    jetPT2darray[2][9][1] = MuonwjetPTresponse29;

    jetPT2darray[3][0][1] = MuonwjetPTresponse30;
    jetPT2darray[3][1][1] = MuonwjetPTresponse31;
    jetPT2darray[3][2][1] = MuonwjetPTresponse32;
    jetPT2darray[3][3][1] = MuonwjetPTresponse33;
    jetPT2darray[3][4][1] = MuonwjetPTresponse34;
    jetPT2darray[3][5][1] = MuonwjetPTresponse35;
    jetPT2darray[3][6][1] = MuonwjetPTresponse36;
    jetPT2darray[3][7][1] = MuonwjetPTresponse37;
    jetPT2darray[3][8][1] = MuonwjetPTresponse38;
    jetPT2darray[3][9][1] = MuonwjetPTresponse39;

    jetPT2darray[4][0][1] = MuonwjetPTresponse40;
    jetPT2darray[4][1][1] = MuonwjetPTresponse41;
    jetPT2darray[4][2][1] = MuonwjetPTresponse42;
    jetPT2darray[4][3][1] = MuonwjetPTresponse43;
    jetPT2darray[4][4][1] = MuonwjetPTresponse44;
    jetPT2darray[4][5][1] = MuonwjetPTresponse45;
    jetPT2darray[4][6][1] = MuonwjetPTresponse46;
    jetPT2darray[4][7][1] = MuonwjetPTresponse47;
    jetPT2darray[4][8][1] = MuonwjetPTresponse48;
    jetPT2darray[4][9][1] = MuonwjetPTresponse49;

    jetPT2darray[5][0][1] = MuonwjetPTresponse50;
    jetPT2darray[5][1][1] = MuonwjetPTresponse51;
    jetPT2darray[5][2][1] = MuonwjetPTresponse52;
    jetPT2darray[5][3][1] = MuonwjetPTresponse53;
    jetPT2darray[5][4][1] = MuonwjetPTresponse54;
    jetPT2darray[5][5][1] = MuonwjetPTresponse55;
    jetPT2darray[5][6][1] = MuonwjetPTresponse56;
    jetPT2darray[5][7][1] = MuonwjetPTresponse57;
    jetPT2darray[5][8][1] = MuonwjetPTresponse58;
    jetPT2darray[5][9][1] = MuonwjetPTresponse59;

    jetPT2darray[6][0][1] = MuonwjetPTresponse60;
    jetPT2darray[6][1][1] = MuonwjetPTresponse61;
    jetPT2darray[6][2][1] = MuonwjetPTresponse62;
    jetPT2darray[6][3][1] = MuonwjetPTresponse63;
    jetPT2darray[6][4][1] = MuonwjetPTresponse64;
    jetPT2darray[6][5][1] = MuonwjetPTresponse65;
    jetPT2darray[6][6][1] = MuonwjetPTresponse66;
    jetPT2darray[6][7][1] = MuonwjetPTresponse67;
    jetPT2darray[6][8][1] = MuonwjetPTresponse68;
    jetPT2darray[6][9][1] = MuonwjetPTresponse69;

    jetPT2darray[7][0][1] = MuonwjetPTresponse70;
    jetPT2darray[7][1][1] = MuonwjetPTresponse71;
    jetPT2darray[7][2][1] = MuonwjetPTresponse72;
    jetPT2darray[7][3][1] = MuonwjetPTresponse73;
    jetPT2darray[7][4][1] = MuonwjetPTresponse74;
    jetPT2darray[7][5][1] = MuonwjetPTresponse75;
    jetPT2darray[7][6][1] = MuonwjetPTresponse76;
    jetPT2darray[7][7][1] = MuonwjetPTresponse77;
    jetPT2darray[7][8][1] = MuonwjetPTresponse78;
    jetPT2darray[7][9][1] = MuonwjetPTresponse79;

    jetPT2darray[8][0][1] = MuonwjetPTresponse80;
    jetPT2darray[8][1][1] = MuonwjetPTresponse81;
    jetPT2darray[8][2][1] = MuonwjetPTresponse82;
    jetPT2darray[8][3][1] = MuonwjetPTresponse83;
    jetPT2darray[8][4][1] = MuonwjetPTresponse84;
    jetPT2darray[8][5][1] = MuonwjetPTresponse85;
    jetPT2darray[8][6][1] = MuonwjetPTresponse86;
    jetPT2darray[8][7][1] = MuonwjetPTresponse87;
    jetPT2darray[8][8][1] = MuonwjetPTresponse88;
    jetPT2darray[8][9][1] = MuonwjetPTresponse89;

    jetPT2darray[9][0][1] = MuonwjetPTresponse90;
    jetPT2darray[9][1][1] = MuonwjetPTresponse91;
    jetPT2darray[9][2][1] = MuonwjetPTresponse92;
    jetPT2darray[9][3][1] = MuonwjetPTresponse93;
    jetPT2darray[9][4][1] = MuonwjetPTresponse94;
    jetPT2darray[9][5][1] = MuonwjetPTresponse95;
    jetPT2darray[9][6][1] = MuonwjetPTresponse96;
    jetPT2darray[9][7][1] = MuonwjetPTresponse97;
    jetPT2darray[9][8][1] = MuonwjetPTresponse98;
    jetPT2darray[9][9][1] = MuonwjetPTresponse99;

    jetPT2darray[0][0][2] = MuonwojetPTresponse00;
    jetPT2darray[0][1][2] = MuonwojetPTresponse01;
    jetPT2darray[0][2][2] = MuonwojetPTresponse02;
    jetPT2darray[0][3][2] = MuonwojetPTresponse03;
    jetPT2darray[0][4][2] = MuonwojetPTresponse04;
    jetPT2darray[0][5][2] = MuonwojetPTresponse05;
    jetPT2darray[0][6][2] = MuonwojetPTresponse06;
    jetPT2darray[0][7][2] = MuonwojetPTresponse07;
    jetPT2darray[0][8][2] = MuonwojetPTresponse08;
    jetPT2darray[0][9][2] = MuonwojetPTresponse09;

    jetPT2darray[1][0][2] = MuonwojetPTresponse10;
    jetPT2darray[1][1][2] = MuonwojetPTresponse11;
    jetPT2darray[1][2][2] = MuonwojetPTresponse12;
    jetPT2darray[1][3][2] = MuonwojetPTresponse13;
    jetPT2darray[1][4][2] = MuonwojetPTresponse14;
    jetPT2darray[1][5][2] = MuonwojetPTresponse15;
    jetPT2darray[1][6][2] = MuonwojetPTresponse16;
    jetPT2darray[1][7][2] = MuonwojetPTresponse17;
    jetPT2darray[1][8][2] = MuonwojetPTresponse18;
    jetPT2darray[1][9][2] = MuonwojetPTresponse19;

    jetPT2darray[2][0][2] = MuonwojetPTresponse20;
    jetPT2darray[2][1][2] = MuonwojetPTresponse21;
    jetPT2darray[2][2][2] = MuonwojetPTresponse22;
    jetPT2darray[2][3][2] = MuonwojetPTresponse23;
    jetPT2darray[2][4][2] = MuonwojetPTresponse24;
    jetPT2darray[2][5][2] = MuonwojetPTresponse25;
    jetPT2darray[2][6][2] = MuonwojetPTresponse26;
    jetPT2darray[2][7][2] = MuonwojetPTresponse27;
    jetPT2darray[2][8][2] = MuonwojetPTresponse28;
    jetPT2darray[2][9][2] = MuonwojetPTresponse29;

    jetPT2darray[3][0][2] = MuonwojetPTresponse30;
    jetPT2darray[3][1][2] = MuonwojetPTresponse31;
    jetPT2darray[3][2][2] = MuonwojetPTresponse32;
    jetPT2darray[3][3][2] = MuonwojetPTresponse33;
    jetPT2darray[3][4][2] = MuonwojetPTresponse34;
    jetPT2darray[3][5][2] = MuonwojetPTresponse35;
    jetPT2darray[3][6][2] = MuonwojetPTresponse36;
    jetPT2darray[3][7][2] = MuonwojetPTresponse37;
    jetPT2darray[3][8][2] = MuonwojetPTresponse38;
    jetPT2darray[3][9][2] = MuonwojetPTresponse39;

    jetPT2darray[4][0][2] = MuonwojetPTresponse40;
    jetPT2darray[4][1][2] = MuonwojetPTresponse41;
    jetPT2darray[4][2][2] = MuonwojetPTresponse42;
    jetPT2darray[4][3][2] = MuonwojetPTresponse43;
    jetPT2darray[4][4][2] = MuonwojetPTresponse44;
    jetPT2darray[4][5][2] = MuonwojetPTresponse45;
    jetPT2darray[4][6][2] = MuonwojetPTresponse46;
    jetPT2darray[4][7][2] = MuonwojetPTresponse47;
    jetPT2darray[4][8][2] = MuonwojetPTresponse48;
    jetPT2darray[4][9][2] = MuonwojetPTresponse49;

    jetPT2darray[5][0][2] = MuonwojetPTresponse50;
    jetPT2darray[5][1][2] = MuonwojetPTresponse51;
    jetPT2darray[5][2][2] = MuonwojetPTresponse52;
    jetPT2darray[5][3][2] = MuonwojetPTresponse53;
    jetPT2darray[5][4][2] = MuonwojetPTresponse54;
    jetPT2darray[5][5][2] = MuonwojetPTresponse55;
    jetPT2darray[5][6][2] = MuonwojetPTresponse56;
    jetPT2darray[5][7][2] = MuonwojetPTresponse57;
    jetPT2darray[5][8][2] = MuonwojetPTresponse58;
    jetPT2darray[5][9][2] = MuonwojetPTresponse59;

    jetPT2darray[6][0][2] = MuonwojetPTresponse60;
    jetPT2darray[6][1][2] = MuonwojetPTresponse61;
    jetPT2darray[6][2][2] = MuonwojetPTresponse62;
    jetPT2darray[6][3][2] = MuonwojetPTresponse63;
    jetPT2darray[6][4][2] = MuonwojetPTresponse64;
    jetPT2darray[6][5][2] = MuonwojetPTresponse65;
    jetPT2darray[6][6][2] = MuonwojetPTresponse66;
    jetPT2darray[6][7][2] = MuonwojetPTresponse67;
    jetPT2darray[6][8][2] = MuonwojetPTresponse68;
    jetPT2darray[6][9][2] = MuonwojetPTresponse69;

    jetPT2darray[7][0][2] = MuonwojetPTresponse70;
    jetPT2darray[7][1][2] = MuonwojetPTresponse71;
    jetPT2darray[7][2][2] = MuonwojetPTresponse72;
    jetPT2darray[7][3][2] = MuonwojetPTresponse73;
    jetPT2darray[7][4][2] = MuonwojetPTresponse74;
    jetPT2darray[7][5][2] = MuonwojetPTresponse75;
    jetPT2darray[7][6][2] = MuonwojetPTresponse76;
    jetPT2darray[7][7][2] = MuonwojetPTresponse77;
    jetPT2darray[7][8][2] = MuonwojetPTresponse78;
    jetPT2darray[7][9][2] = MuonwojetPTresponse79;

    jetPT2darray[8][0][2] = MuonwojetPTresponse80;
    jetPT2darray[8][1][2] = MuonwojetPTresponse81;
    jetPT2darray[8][2][2] = MuonwojetPTresponse82;
    jetPT2darray[8][3][2] = MuonwojetPTresponse83;
    jetPT2darray[8][4][2] = MuonwojetPTresponse84;
    jetPT2darray[8][5][2] = MuonwojetPTresponse85;
    jetPT2darray[8][6][2] = MuonwojetPTresponse86;
    jetPT2darray[8][7][2] = MuonwojetPTresponse87;
    jetPT2darray[8][8][2] = MuonwojetPTresponse88;
    jetPT2darray[8][9][2] = MuonwojetPTresponse89;

    jetPT2darray[9][0][2] = MuonwojetPTresponse90;
    jetPT2darray[9][1][2] = MuonwojetPTresponse91;
    jetPT2darray[9][2][2] = MuonwojetPTresponse92;
    jetPT2darray[9][3][2] = MuonwojetPTresponse93;
    jetPT2darray[9][4][2] = MuonwojetPTresponse94;
    jetPT2darray[9][5][2] = MuonwojetPTresponse95;
    jetPT2darray[9][6][2] = MuonwojetPTresponse96;
    jetPT2darray[9][7][2] = MuonwojetPTresponse97;
    jetPT2darray[9][8][2] = MuonwojetPTresponse98;
    jetPT2darray[9][9][2] = MuonwojetPTresponse99;

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
                    jetPTresponse -> Fill(JERtmp);
                    /*
		    for (Int_t muonentry=0;muonentry < nMuon; muonentry++) {
                    	recoMuonEta = MuonEta->GetValue(muonentry);
                    	recoMuonPhi = MuonPhi->GetValue(muonentry);
                    	recoMuonPt = MuonPt->GetValue(muonentry);
			recoMuon.SetPtEtaPhiM(recoMuonPt, recoMuonEta, recoMuonPhi, 0);
			MuonDeltaPhi = recoMuon.DeltaPhi(AKTjet);
			MuonDeltaR = recoMuon.DeltaR(AKTjet);
		   	if (MuonDeltaR < 0.5) {
		      	    MuonTagging = true;
			    //AKTjet = AKTjet + recoMuon;
			    JERtmp = AKTjet.Pt()/GenJet.Pt();
		            MuonwjetPTresponse->Fill(JERtmp);
			    break;
		        }
		    }
		    */

                    /*
		    for (Int_t partentry = 0; partentry < nParticle; partentry++) {
                        ParticlePID = PID -> GetValue(partentry);
                        if (TMath::Abs(ParticlePID) == 13) {
                            GenMuonEta = GenParticleEta -> GetValue(partentry);
                            GenMuonPhi = GenParticlePhi -> GetValue(partentry);
                            GenMuonPt = GenParticlePt -> GetValue(partentry);
                            GenMuon.SetPtEtaPhiM(GenMuonPt, GenMuonEta, GenMuonPhi, 0);

                            GenMuonDeltaPhi = TMath::Abs(GenMuonPhi - MatchedGenJetPhi);
                            if (GenMuonDeltaPhi > TMath::Pi()) {
                                GenMuonDeltaPhi -= TMath::TwoPi();
                            }
                            GenMuonDeltaR = TMath::Sqrt(pow((GenMuonEta - MatchedGenJetEta), 2) + pow(GenMuonDeltaPhi, 2));
                            //GenMuonDeltaR = GenMuon.DeltaR(AKTjet);
                            if (GenMuonDeltaR < 0.5) {
                                MuonTagging = true;
                                //AKTjet += GenMuon;
                                JERtmp = AKTjet.Pt() / GenJet.Pt();
                                break;
                            }
                        }
                    }
		    */

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
		    /*
                    if (MuonTagging == true) {
                        matchJetMuonwCnt[ThetaGrid][PTGrid] += 1;
                        MuonwJER[ThetaGrid][PTGrid] = (MuonwJER[ThetaGrid][PTGrid] *(matchJetMuonwCnt[ThetaGrid][PTGrid] - 1) + JERtmp) / matchJetMuonwCnt[ThetaGrid][PTGrid];
                        jetPT2darray[ThetaGrid][PTGrid][1] -> Fill(JERtmp);
                        MuonwjetPTresponse -> Fill(JERtmp);

                    } else {
                        matchJetMuonwoCnt[ThetaGrid][PTGrid] += 1;
                        MuonwoJER[ThetaGrid][PTGrid] = (MuonwoJER[ThetaGrid][PTGrid] *(matchJetMuonwoCnt[ThetaGrid][PTGrid] - 1) + JERtmp) / matchJetMuonwoCnt[ThetaGrid][PTGrid];
                        jetPT2darray[ThetaGrid][PTGrid][2] -> Fill(JERtmp);
                    }
		    */
                    jetPT2darray[ThetaGrid][PTGrid][0] -> Fill(JERtmp);
                } else {
                    //cout << DeltaR << endl;
                }
            }
        }
    }
    cout << endl << "Calibration's result:" << endl;
    //Calibration result
    TCanvas *myCaliCanvas = new TCanvas("myCaliCanvas", "My Canvas", 800, 600);
    jetPTresponse -> Draw();
    //myCaliCanvas -> SaveAs("jetPTresponse.png");
    
    /*
    MuonwjetPTresponse -> Draw();
    myCaliCanvas -> SaveAs("MuonwjetPTresponse.png");
    */
    for (Int_t PTentry = 9; PTentry >= 0; PTentry -= 1) {
        for (Int_t ThetaEntry = 0; ThetaEntry < 10; ThetaEntry++) {
            if (matchJetCnt[ThetaEntry][PTentry] == 0) {
                JER[ThetaEntry][PTentry] = 1.0000;
            }
	    /*
            if (matchJetMuonwCnt[ThetaEntry][PTentry] == 0) {
                MuonwJER[ThetaEntry][PTentry] = 1.0000;
            }
            if (matchJetMuonwoCnt[ThetaEntry][PTentry] == 0) {
                MuonwoJER[ThetaEntry][PTentry] = 1.0000;
            }
	    */
            jetPT2darray[ThetaEntry][PTentry][0] -> GetXaxis() -> SetTitle("Jet P_{T} Response");
            jetPT2darray[ThetaEntry][PTentry][0] -> Draw();
            myCaliCanvas -> SaveAs(TString::Format("jetPTresponse_%d%d.png", ThetaEntry, PTentry));
            jetPT2darray[ThetaEntry][PTentry][0] -> Write();
            /*
	    jetPT2darray[ThetaEntry][PTentry][1] -> GetXaxis() -> SetTitle("Jet P_{T} Response");
            jetPT2darray[ThetaEntry][PTentry][1] -> Draw();
            myCaliCanvas -> SaveAs(TString::Format("MuonwjetPTresponse_%d%d.png", ThetaEntry, PTentry));
            jetPT2darray[ThetaEntry][PTentry][1] -> Write();
            jetPT2darray[ThetaEntry][PTentry][2] -> GetXaxis() -> SetTitle("Jet P_{T} Response");
            jetPT2darray[ThetaEntry][PTentry][2] -> Draw();
            myCaliCanvas -> SaveAs(TString::Format("MuonwojetPTresponse_%d%d.png", ThetaEntry, PTentry));
            jetPT2darray[ThetaEntry][PTentry][2] -> Write();
	    */
        }
    }
    cout << "Jet Energy Calibration gives following Jet Energy Response:" << endl;
    for (Int_t PTentry = 9; PTentry >= 0; PTentry -= 1) {
        for (Int_t ThetaEntry = 0; ThetaEntry < 10; ThetaEntry++) {
            cout << floorf(JER[ThetaEntry][PTentry] *100) / 100 << ", ";
        }
        cout << endl;
    }
    /*
    cout << endl << "MuonTagged" << endl;
    for (Int_t PTentry = 9; PTentry >= 0; PTentry -= 1) {
        for (Int_t ThetaEntry = 0; ThetaEntry < 10; ThetaEntry++) {
            cout << floorf(MuonwJER[ThetaEntry][PTentry] *100) / 100 << ", ";
        }
        cout << endl;
    }
    cout << endl << "Not MuonTagged" << endl;
    for (Int_t PTentry = 9; PTentry >= 0; PTentry -= 1) {
        for (Int_t ThetaEntry = 0; ThetaEntry < 10; ThetaEntry++) {
            cout << floorf(MuonwoJER[ThetaEntry][PTentry] *100) / 100 << ", ";
        }
        cout << endl;
    }
    */
    cout << endl << "Number of Jets included in calculation" << endl;
    for (Int_t PTentry = 9; PTentry >= 0; PTentry -= 1) {
        for (Int_t ThetaEntry = 0; ThetaEntry < 10; ThetaEntry++) {
            cout << matchJetCnt[ThetaEntry][PTentry] << ", ";
        }
        cout << endl;
    }
    /*
    cout << endl << "MuonTagged" << endl;
    for (Int_t PTentry = 9; PTentry >= 0; PTentry -= 1) {
        for (Int_t ThetaEntry = 0; ThetaEntry < 10; ThetaEntry++) {
            if (matchJetCnt[ThetaEntry][PTentry] == 0) {
                cout << "N/A" << ", ";
            } else {
                cout << floorf(100 *matchJetMuonwCnt[ThetaEntry][PTentry] / matchJetCnt[ThetaEntry][PTentry]) << "%, ";
            }
            //cout << matchJetMuonwCnt[ThetaEntry][PTentry] << ", ";
        }
        cout << endl;
    }
    cout << endl << "Not MuonTagged" << endl;
    for (Int_t PTentry = 9; PTentry >= 0; PTentry -= 1) {
        for (Int_t ThetaEntry = 0; ThetaEntry < 10; ThetaEntry++) {
            if (matchJetCnt[ThetaEntry][PTentry] == 0) {
                cout << "N/A" << ", ";
            } else {
                cout << floorf(100 *matchJetMuonwoCnt[ThetaEntry][PTentry] / matchJetCnt[ThetaEntry][PTentry]) << "%, ";
            }
            //cout << matchJetMuonwoCnt[ThetaEntry][PTentry] << ", ";
        }
        cout << endl;
    }
    */

    tree_output -> Write();
    output -> Close();
    file_sig -> Close();

    cout << endl << "Calibration all done! Quit..." << endl;
}

void Pairing_w_JES(const char *inputFile,
    Double_t JER[10][10],
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

    TH2D *jet1Reso_Pt = new TH2D("jet1Reso_Pt", "jet1Reso_Pt", 70, 0, 400, 70, -1, 4);
    TH2D *jet2Reso_Pt = new TH2D("jet2Reso_Pt", "jet2Reso_Pt", 70, 0, 400, 70, -1, 4);
    TH2D *jet3Reso_Pt = new TH2D("jet3Reso_Pt", "jet3Reso_Pt", 70, 0, 400, 70, -1, 4);
    TH2D *jet4Reso_Pt = new TH2D("jet4Reso_Pt", "jet4Reso_Pt", 70, 0, 400, 70, -1, 4);

    TH1D *badjet1theta = new TH1D("badjet1theta", "badjet1theta", 50, 0, 3.1415);
    TH1D *badjet2theta = new TH1D("badjet2theta", "badjet2theta", 50, 0, 3.1315);
    TH1D *badjet3theta = new TH1D("badjet3theta", "badjet3theta", 50, 0, 3.1415);
    TH1D *badjet4theta = new TH1D("badjet4theta", "badjet4theta", 50, 0, 3.1415);

    TH1D *goodjet1theta = new TH1D("goodjet1theta", "goodjet1theta", 50, 0, 3.1415);
    TH1D *goodjet2theta = new TH1D("goodjet2theta", "goodjet2theta", 50, 0, 3.1415);
    TH1D *goodjet3theta = new TH1D("goodjet3theta", "goodjet3theta", 50, 0, 3.1415);
    TH1D *goodjet4theta = new TH1D("goodjet4theta", "goodjet4theta", 50, 0, 3.1415);

    TH1D *alljet1theta = new TH1D("alljet1theta", "alljet1theta", 50, 0, 3.1415);
    TH1D *alljet2theta = new TH1D("alljet2theta", "alljet2theta", 50, 0, 3.1315);
    TH1D *alljet3theta = new TH1D("alljet3theta", "alljet3theta", 50, 0, 3.1415);
    TH1D *alljet4theta = new TH1D("alljet4theta", "alljet4theta", 50, 0, 3.1415);

    TH2D *jet1Reso_DeltaR = new TH2D("jet1Reso_DeltaR", "jet1Reso_DeltaR", 50, 0, 0.5, 50, -1, 4);
    TH2D *jet2Reso_DeltaR = new TH2D("jet2Reso_DeltaR", "jet2Reso_DeltaR", 50, 0, 0.5, 50, -1, 4);
    TH2D *jet3Reso_DeltaR = new TH2D("jet3Reso_DeltaR", "jet3Reso_DeltaR", 50, 0, 0.5, 50, -1, 4);
    TH2D *jet4Reso_DeltaR = new TH2D("jet4Reso_DeltaR", "jet4Reso_DeltaR", 50, 0, 0.5, 50, -1, 4);

    TH1D *badjet1DeltaR = new TH1D("badjet1DeltaR", "badjet1DeltaR", 50, 0, 0.5);
    TH1D *badjet2DeltaR = new TH1D("badjet2DeltaR", "badjet2DeltaR", 50, 0, 0.5);
    TH1D *badjet3DeltaR = new TH1D("badjet3DeltaR", "badjet3DeltaR", 50, 0, 0.5);
    TH1D *badjet4DeltaR = new TH1D("badjet4DeltaR", "badjet4DeltaR", 50, 0, 0.5);

    TH1D *alljet1DeltaR = new TH1D("alljet1DeltaR", "alljet1DeltaR", 50, 0, 0.5);
    TH1D *alljet2DeltaR = new TH1D("alljet2DeltaR", "alljet2DeltaR", 50, 0, 0.5);
    TH1D *alljet3DeltaR = new TH1D("alljet3DeltaR", "alljet3DeltaR", 50, 0, 0.5);
    TH1D *alljet4DeltaR = new TH1D("alljet4DeltaR", "alljet4DeltaR", 50, 0, 0.5);

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

    Int_t AKTjet1BTag1;
    Int_t AKTjet1BTag2;
    Int_t AKTjet2BTag1;
    Int_t AKTjet2BTag2;

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
    TLorentzVector AKT1jet1;
    TLorentzVector AKT1jet2;
    TLorentzVector AKT2jet1;
    TLorentzVector AKT2jet2;

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

    Int_t BDTjet1BTag1;
    Int_t BDTjet1BTag2;
    Int_t BDTjet2BTag1;
    Int_t BDTjet2BTag2;

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

    Double_t DiHiggsMETDeltaTheta;
    Int_t NBTag;
    tree_BDT -> Branch("BDTNjets", &BDTNjets, "BDTNjets/I");
    tree_BDT -> Branch("NBTag", &NBTag, "NBTag/I");

    tree_BDT -> Branch("BDTjet1BTag1", &BDTjet1BTag1, "BDTjet1BTag1/I");
    tree_BDT -> Branch("BDTjet1BTag2", &BDTjet1BTag2, "BDTjet1BTag2/I");
    tree_BDT -> Branch("BDTjet2BTag1", &BDTjet2BTag1, "BDTjet2BTag1/I");
    tree_BDT -> Branch("BDTjet2BTag2", &BDTjet2BTag2, "BDTjet2BTag2/I");

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
        tree_sig -> GetEntry(entry);
        //tree_output -> GetEntry(entry);
        AKTjet_size -> GetBranch() -> GetEntry(entry);
        GenJet_size -> GetBranch() -> GetEntry(entry);

        Int_t nAKTjet = AKTjet_size -> GetValue();
        Int_t nGenJet = GenJet_size -> GetValue();

        AKTjet_eta -> GetBranch() -> GetEntry(entry);
        AKTjet_phi -> GetBranch() -> GetEntry(entry);
        AKTjet_pt -> GetBranch() -> GetEntry(entry);
        AKTjet_mass -> GetBranch() -> GetEntry(entry);
        AKTjet_BTag -> GetBranch() -> GetEntry(entry);

        GenJet_eta -> GetBranch() -> GetEntry(entry);
        GenJet_phi -> GetBranch() -> GetEntry(entry);
        GenJet_pt -> GetBranch() -> GetEntry(entry);
        GenJet_mass -> GetBranch() -> GetEntry(entry);

	MissingETMET -> GetBranch() -> GetEntry(entry);
	MissingETEta -> GetBranch() -> GetEntry(entry);
	MissingETPhi -> GetBranch() -> GetEntry(entry);

	BDTMissingETMET = MissingETMET -> GetValue();
	BDTMissingETEta = MissingETEta -> GetValue();
	BDTMissingETPhi = MissingETPhi -> GetValue();

        AKTjetpair1Mass = 10000;
        AKTjetpair2Mass = 10000;

        pair1jet1entry = 0;
        pair1jet2entry = 0;
        pair2jet1entry = 0;
        pair2jet2entry = 0;

        AKTjet2entry1 = 0;
        AKTjet2entry2 = 0;

        diHiggsdis = 10000;

        //Pairing up Anti-kt jets

        Double_t AKTjetpair[nAKTjet *(nAKTjet - 1) / 2][6];

        AKTjetpaircnt = 0;

        if (nAKTjet >= 4) {
            for (Int_t aktentry = 0; aktentry < nAKTjet; aktentry++) {
                AKTjetEta = AKTjet_eta -> GetValue(aktentry);
                AKTjetPt = AKTjet_pt -> GetValue(aktentry);
                AKTjetTheta = 2 *atan(exp(-AKTjetEta));
                AKTjetPT_Theta -> Fill(AKTjetTheta, AKTjetPt);
                JetEnergyFix(AKTjetEta, AKTjetPt, JER);
            }
            for (Int_t genentry = 0; genentry < nGenJet; genentry++) {
                GenJetEta = GenJet_eta -> GetValue(genentry);
                GenJetPt = GenJet_pt -> GetValue(genentry);
                GenJetTheta = 2 *atan(exp(-GenJetEta));
                GenJetPT_Theta -> Fill(GenJetTheta, GenJetPt);
            }
            //Pairing up leading anti-KT jet pair
            for (Int_t akt1entry = 0; akt1entry < nAKTjet; akt1entry++) {
                for (Int_t akt2entry = akt1entry + 1; akt2entry < nAKTjet; akt2entry++) {

                    AKTjet1eta1 = AKTjet_eta -> GetValue(akt1entry);
                    AKTjet1phi1 = AKTjet_phi -> GetValue(akt1entry);
                    AKTjet1pt1 = AKTjet_pt -> GetValue(akt1entry);
                    JetEnergyFix(AKTjet1eta1, AKTjet1pt1, JER);
                    AKTjet1mass1 = AKTjet_mass -> GetValue(akt1entry);
                    AKTjet1eta2 = AKTjet_eta -> GetValue(akt2entry);
                    AKTjet1phi2 = AKTjet_phi -> GetValue(akt2entry);
                    AKTjet1pt2 = AKTjet_pt -> GetValue(akt2entry);
                    JetEnergyFix(AKTjet1eta2, AKTjet1pt2, JER);
                    AKTjet1mass2 = AKTjet_mass -> GetValue(akt2entry);
                    AKTjetpairmass = 0;

                    AKT1jet1.SetPtEtaPhiM(AKTjet1pt1, AKTjet1eta1, AKTjet1phi1, AKTjet1mass1);
                    AKT1jet2.SetPtEtaPhiM(AKTjet1pt2, AKTjet1eta2, AKTjet1phi2, AKTjet1mass2);
                    AKTh1 = AKT1jet1 + AKT1jet2;
                    AKTjetpairmass = AKTh1.Mag();

                    AKTjetpair[AKTjetpaircnt][0] = AKTjetpairmass;
                    AKTjetpair[AKTjetpaircnt][1] = akt1entry;
                    AKTjetpair[AKTjetpaircnt][2] = akt2entry;

                    Int_t akt1entry2;
                    Int_t akt2entry2;
                    for (akt1entry2 = 0; akt1entry2 < nAKTjet; akt1entry2++) {
                        if ((akt1entry2 != akt1entry) and(akt1entry2 != akt2entry)) {
                            for (akt2entry2 = 0; akt2entry2 < nAKTjet; akt2entry2++) {
                                if (((akt2entry2 != akt1entry) and(akt2entry2 != akt2entry)) and(akt2entry2 != akt1entry)) {

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
                                    AKTjetpairmass = 0;

                                    AKT2jet1.SetPtEtaPhiM(AKTjet2pt1, AKTjet2eta1, AKTjet2phi1, AKTjet2mass1);
                                    AKT2jet2.SetPtEtaPhiM(AKTjet2pt2, AKTjet2eta2, AKTjet2phi2, AKTjet2mass2);
                                    AKTh2 = AKT2jet1 + AKT2jet2;
                                    AKTjetpairmass = AKTh2.Mag();

                                    if (abs(125 - AKTjetpairmass) < abs(125 - AKTjetpair2Mass)) {
                                        AKTjetpair2Mass = AKTjetpairmass;
                                        AKTjet2entry1 = akt1entry2;
                                        AKTjet2entry2 = akt2entry2;
                                    }
                                }
                            }
                        }
                    }
                    AKTjetpair[AKTjetpaircnt][3] = AKTjetpair2Mass;
                    AKTjetpair[AKTjetpaircnt][4] = AKTjet2entry1;
                    AKTjetpair[AKTjetpaircnt][5] = AKTjet2entry2;
                }

                AKTjetpaircnt += 1;
            }
            //cout << AKTjetpaircnt;
            for (Int_t entry = 0; entry < (nAKTjet - 1) *nAKTjet / 2; entry++) {
                if (pow((125 - AKTjetpair[entry][0]), 2) + pow((125 - AKTjetpair[entry][3]), 2) < diHiggsdis) {

                    pair1jet1entry = AKTjetpair[entry][1];
                    pair1jet2entry = AKTjetpair[entry][2];
                    pair2jet1entry = AKTjetpair[entry][4];
                    pair2jet2entry = AKTjetpair[entry][5];

                    AKTjetpair1Mass = AKTjetpair[entry][0];
                    AKTjetpair2Mass = AKTjetpair[entry][3];

                    //diHiggsdis = abs(125 - AKTjetpair[entry][0]) *abs(125 - AKTjetpair[entry][3]);
                    diHiggsdis = pow((125 - AKTjetpair[entry][0]), 2) + pow((125 - AKTjetpair[entry][3]), 2);
                }
            }

            AKTjet1eta1 = AKTjet_eta -> GetValue(pair1jet1entry);
            AKTjet1phi1 = AKTjet_phi -> GetValue(pair1jet1entry);
            AKTjet1pt1 = AKTjet_pt -> GetValue(pair1jet1entry);
            JetEnergyFix(AKTjet1eta1, AKTjet1pt1, JER);
            AKTjet1mass1 = AKTjet_mass -> GetValue(pair1jet1entry);
            AKTjet1eta2 = AKTjet_eta -> GetValue(pair1jet2entry);
            AKTjet1phi2 = AKTjet_phi -> GetValue(pair1jet2entry);
            AKTjet1pt2 = AKTjet_pt -> GetValue(pair1jet2entry);
            JetEnergyFix(AKTjet1eta2, AKTjet1pt2, JER);
            AKTjet1mass2 = AKTjet_mass -> GetValue(pair1jet2entry);

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

	    AKTjet1BTag1 = AKTjet_BTag -> GetValue(pair1jet1entry);
	    AKTjet1BTag2 = AKTjet_BTag -> GetValue(pair1jet2entry);
	    AKTjet2BTag1 = AKTjet_BTag -> GetValue(pair2jet1entry);
	    AKTjet2BTag2 = AKTjet_BTag -> GetValue(pair2jet2entry);

            bool AKT1jet1flag = false;
            bool AKT1jet2flag = false;
            bool AKT2jet1flag = false;
            bool AKT2jet2flag = false;

            bool AKT1jet1edgeflag = false;
            bool AKT1jet2edgeflag = false;
            bool AKT2jet1edgeflag = false;
            bool AKT2jet2edgeflag = false;

            bool disCutFlag = false;

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
            //Truth matching for anti-KT jets
            for (Int_t gen1entry = 0; gen1entry < nGenJet; gen1entry++) {
                Gen1eta = GenJet_eta -> GetValue(gen1entry);
                Gen1phi = GenJet_phi -> GetValue(gen1entry);
                jet1DeltaPhi = std::abs(AKTjet1phi1 - Gen1phi);
                if (jet1DeltaPhi > TMath::Pi()) {
                    jet1DeltaPhi -= 2 *TMath::Pi();
                }
                Double_t jet1DeltaRtmp = TMath::Sqrt(pow((AKTjet1eta1 - Gen1eta), 2) + pow(jet1DeltaPhi, 2));
                if (jet1DeltaRtmp < jet1DeltaR) {
                    jet1DeltaR = jet1DeltaRtmp;
                    if (jet1DeltaR < 0.5) {
                        //if (abs(Gen1eta) < 2.25){
                        AKT1jet1flag = true;
                        jet1entry = gen1entry;
                        //}
                    }
                }
            }
            for (Int_t gen2entry = 0; gen2entry < nGenJet; gen2entry++) {
                Gen2eta = GenJet_eta -> GetValue(gen2entry);
                Gen2phi = GenJet_phi -> GetValue(gen2entry);
                jet2DeltaPhi = std::abs(AKTjet1phi2 - Gen2phi);
                if (jet2DeltaPhi > TMath::Pi()) {
                    jet2DeltaPhi -= 2 *TMath::Pi();
                }
                Double_t jet2DeltaRtmp = TMath::Sqrt(pow((AKTjet1eta2 - Gen2eta), 2) + pow(jet2DeltaPhi, 2));
                if (gen2entry != jet1entry) {
                    if (jet2DeltaRtmp < jet2DeltaR) {
                        jet2DeltaR = jet2DeltaRtmp;
                        if (jet2DeltaR < 0.5) {
                            //if (abs(Gen2eta) < 2.25){
                            AKT1jet2flag = true;
                            jet2entry = gen2entry;
                            //}
                        }
                    }
                }
            }

            for (Int_t gen3entry = 0; gen3entry < nGenJet; gen3entry++) {
                Gen3eta = GenJet_eta -> GetValue(gen3entry);
                Gen3phi = GenJet_phi -> GetValue(gen3entry);
                jet3DeltaPhi = std::abs(AKTjet2phi1 - Gen3phi);
                if (jet3DeltaPhi > TMath::Pi()) {
                    jet3DeltaPhi -= 2 *TMath::Pi();
                }
                Double_t jet1DeltaRtmp = TMath::Sqrt(pow((AKTjet2eta1 - Gen3eta), 2) + pow(jet3DeltaPhi, 2));
                if (jet1DeltaRtmp < jet3DeltaR) {
                    jet3DeltaR = jet1DeltaRtmp;
                    if (jet3DeltaR < 0.5) {
                        //if (abs(Gen1eta) < 2.25){
                        AKT2jet1flag = true;
                        jet3entry = gen3entry;
                        //}
                    }
                }
            }

            for (Int_t gen4entry = 0; gen4entry < nGenJet; gen4entry++) {
                Gen4eta = GenJet_eta -> GetValue(gen4entry);
                Gen4phi = GenJet_phi -> GetValue(gen4entry);
                jet4DeltaPhi = std::abs(AKTjet2phi2 - Gen4phi);
                if (jet4DeltaPhi > TMath::Pi()) {
                    jet4DeltaPhi -= 2 *TMath::Pi();
                }
                Double_t jet2DeltaRtmp = TMath::Sqrt(pow((AKTjet2eta2 - Gen4eta), 2) + pow(jet4DeltaPhi, 2));
                if (gen4entry != jet1entry) {
                    if (jet2DeltaRtmp < jet4DeltaR) {
                        jet4DeltaR = jet2DeltaRtmp;
                        if (jet4DeltaR < 0.5) {
                            //if (abs(Gen2eta) < 2.25){
                            AKT2jet2flag = true;
                            jet4entry = gen4entry;
                            //}
                        }
                    }
                }
            }


	    if (/*TMath::Abs(AKTjetpair1Mass - 125) < 20 and TMath::Abs(AKTjetpair2Mass - 125) < 20 and*/ TMath::Abs(AKTjetpair1Mass-AKTjetpair2Mass) < 20) {
		disCutFlag = true;
	    }

            if (AKTjetpair1Mass < AKTjetpair2Mass) {
                swap(AKTjetpair1Mass, AKTjetpair2Mass);
                swap(jet1entry, jet3entry);
                swap(jet2entry, jet4entry);
                swap(AKTjet1pt1, AKTjet2pt1);
                swap(AKTjet1pt2, AKTjet2pt2);
                swap(AKTjet1eta1, AKTjet2eta1);
		swap(AKTjet1eta2, AKTjet2eta2);
                swap(AKTjet1phi1, AKTjet2phi1);
                swap(AKTjet1phi2, AKTjet2phi2);
		swap(AKTjet1mass1, AKTjet2mass1);
		swap(AKTjet1mass2, AKTjet2mass2);
		swap(AKTjet1BTag1, AKTjet2BTag1);
		swap(AKTjet1BTag2, AKTjet2BTag2);
            }

            Gen1eta = GenJet_eta -> GetValue(jet1entry);
            Gen1phi = GenJet_phi -> GetValue(jet1entry);
            Gen1pt = GenJet_pt -> GetValue(jet1entry);
            Gen1mass = GenJet_mass -> GetValue(jet1entry);
            Gen2eta = GenJet_eta -> GetValue(jet2entry);
            Gen2phi = GenJet_phi -> GetValue(jet2entry);
            Gen2pt = GenJet_pt -> GetValue(jet2entry);
            Gen2mass = GenJet_mass -> GetValue(jet2entry);

            GenJet1.SetPtEtaPhiM(Gen1pt, Gen1eta, Gen1phi, Gen1mass);
            GenJet2.SetPtEtaPhiM(Gen2pt, Gen2eta, Gen2phi, Gen2mass);
            Genh2 = GenJet1 + GenJet2;
            GenJetMass = Genh2.Mag();
            Double_t AKTGenMass1diff = (AKTjetpair1Mass - GenJetMass) / GenJetMass;

            Gen3eta = GenJet_eta -> GetValue(jet3entry);
            Gen3phi = GenJet_phi -> GetValue(jet3entry);
            Gen3pt = GenJet_pt -> GetValue(jet3entry);
            Gen3mass = GenJet_mass -> GetValue(jet3entry);
            Gen4eta = GenJet_eta -> GetValue(jet4entry);
            Gen4phi = GenJet_phi -> GetValue(jet4entry);
            Gen4pt = GenJet_pt -> GetValue(jet4entry);
            Gen4mass = GenJet_mass -> GetValue(jet4entry);

            GenJet3.SetPtEtaPhiM(Gen3pt, Gen3eta, Gen3phi, Gen3mass);
            GenJet4.SetPtEtaPhiM(Gen4pt, Gen4eta, Gen4phi, Gen4mass);
            Genh2 = GenJet3 + GenJet4;
            GenJetMass = Genh2.Mag();
            Double_t AKTGenMass2diff = (AKTjetpair2Mass - GenJetMass) / GenJetMass;

            Double_t AKTpt[] = {
                AKTjet1pt1,
                AKTjet1pt2,
                AKTjet2pt1,
                AKTjet2pt2
            };
            //int AKTsize = sizeof(AKTpt) / sizeof(AKTpt[0]);

            Double_t Genpt[] = {
                Gen1pt,
                Gen2pt,
                Gen3pt,
                Gen4pt
            };
            //int Gensize = sizeof(Genpt) / sizeof(Genpt[0]);

            Double_t deltaR[] = {
                jet1DeltaR,
                jet2DeltaR,
                jet3DeltaR,
                jet4DeltaR
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


			    if (AKTBTag[0] >=4){
				BDTjet1BTag1 = 1;
			    } else {
				BDTjet1BTag1 = 0;
			    }

			    if (AKTBTag[1] >=4){
				BDTjet1BTag2 = 1;
			    } else {
				BDTjet1BTag2 = 0;
			    }

			    if (AKTBTag[2] >=4){
				BDTjet2BTag1 = 1;
			    } else {
				BDTjet2BTag1 = 0;
			    }

			    if (AKTBTag[3] >=4){
				BDTjet2BTag2 = 1;
			    } else {
				BDTjet2BTag2 = 0;
			    }

			    NBTag = BDTjet1BTag1 + BDTjet1BTag2 + BDTjet2BTag1 + BDTjet2BTag2;

            /*
            //sort pt
            for (int i = 0; i < 3; i++) {
                for (int j = i + 1; j < 4; j++) {
                    if (AKTpt[j] > AKTpt[i]) {
                        swap(AKTpt[j], AKTpt[i]);
                        swap(Genpt[j], Genpt[i]);
                        swap(deltaR[j], deltaR[i]);
                        swap(AKTeta[j], AKTeta[i]);
			swap(AKTphi[j], AKTphi[i]);
			swap(AKTmass[j], AKTmass[i]);
                    }
                }
            }
	    */
            //calculate resolution
            Double_t AKTGenPt1diff = (AKTpt[0] - Genpt[0]) / Genpt[0];
            Double_t AKTGenPt2diff = (AKTpt[1] - Genpt[1]) / Genpt[1];
            Double_t AKTGenPt3diff = (AKTpt[2] - Genpt[2]) / Genpt[2];
            Double_t AKTGenPt4diff = (AKTpt[3] - Genpt[3]) / Genpt[3];

            //calculate angle
            AKTjet1theta1 = 2 *atan(exp(-AKTeta[0]));
            AKTjet1theta2 = 2 *atan(exp(-AKTeta[1]));
            AKTjet2theta1 = 2 *atan(exp(-AKTeta[2]));
            AKTjet2theta2 = 2 *atan(exp(-AKTeta[3]));

            //theta edge check
            /*
	     if (AKTjet1theta1>0.5 and AKTjet1theta1<2.5){
	         AKT1jet1edgeflag = true;
	     } else {
	         AKT1jet1edgeflag = false;
	     }

             if (AKTjet1theta2>0.5 and AKTjet1theta2<2.5){
	         AKT1jet2edgeflag = true;
	     } else {
	         AKT1jet2edgeflag = false;
	     }

	     if (AKTjet2theta1>0.5 and AKTjet2theta1<2.5){
	         AKT2jet1edgeflag = true;
	     } else {
	         AKT2jet1edgeflag = false;
	     }

	     if (AKTjet2theta2>0.5 and AKTjet2theta2<2.5){
	         AKT2jet2edgeflag = true;
	     } else {
	         AKT2jet2edgeflag = false;
	     }
         */

            //if (AKT1jet1flag == true and AKT1jet2flag == true and AKT2jet1flag == true and AKT2jet2flag == true){
	    if (NBTag>=4){
            	AKTGenPt1Comp -> Fill(AKTGenPt1diff);
            	AKTGenPt2Comp -> Fill(AKTGenPt2diff);
            	AKTGenPt3Comp -> Fill(AKTGenPt3diff);
            	AKTGenPt4Comp -> Fill(AKTGenPt4diff);
            	alljet1DeltaR -> Fill(deltaR[0]);
            	alljet2DeltaR -> Fill(deltaR[1]);
            	alljet3DeltaR -> Fill(deltaR[2]);
            	alljet4DeltaR -> Fill(deltaR[3]);
            	alljet1theta -> Fill(AKTjet1theta1);
            	alljet2theta -> Fill(AKTjet1theta2);
            	alljet3theta -> Fill(AKTjet2theta1);
            	alljet4theta -> Fill(AKTjet2theta2);
            	if (abs(AKTGenPt1diff) > 0.2) {
                    jet1Reso_Pt -> Fill(Genpt[0], AKTGenPt1diff);
                    jet1Reso_DeltaR -> Fill(deltaR[0], AKTGenPt1diff);
                    badjet1DeltaR -> Fill(deltaR[0]);
                    badjet1theta -> Fill(AKTjet1theta1);
            	} else {
                    goodjet1theta -> Fill(AKTjet1theta1);
            	}
            	if (abs(AKTGenPt2diff) > 0.2) {
                    jet2Reso_Pt -> Fill(Genpt[1], AKTGenPt2diff);
                    jet2Reso_DeltaR -> Fill(deltaR[1], AKTGenPt2diff);
                    badjet2DeltaR -> Fill(deltaR[1]);
                    badjet2theta -> Fill(AKTjet1theta2);
            	} else {
                    goodjet2theta -> Fill(AKTjet1theta2);
            	}
            	if (abs(AKTGenPt3diff) > 0.2) {
                    jet3Reso_Pt -> Fill(Genpt[2], AKTGenPt3diff);
                    jet3Reso_DeltaR -> Fill(deltaR[2], AKTGenPt3diff);
                    badjet3DeltaR -> Fill(deltaR[2]);
                    badjet3theta -> Fill(AKTjet2theta1);
            	} else {
                    goodjet3theta -> Fill(AKTjet2theta1);
            	}
            	if (abs(AKTGenPt4diff) > 0.2) {
                    jet4Reso_Pt -> Fill(Genpt[3], AKTGenPt4diff);
                    jet4Reso_DeltaR -> Fill(deltaR[3], AKTGenPt4diff);
                    badjet4DeltaR -> Fill(deltaR[3]);
                    badjet4theta -> Fill(AKTjet2theta2);
            	} else {
                    goodjet4theta -> Fill(AKTjet2theta2);
            	}

            	//if (abs(AKTGenPt1diff) < 0.15 and abs(AKTGenPt2diff) < 0.15 and abs(AKTGenPt3diff) < 0.15 and abs(AKTGenPt4diff) < 0.15) {
                    //if (AKT1jet1edgeflag==true and AKT1jet2edgeflag==true and AKT2jet1edgeflag==true and AKT2jet2edgeflag==true) {
		        //if (disCutFlag == true) {
	   	  	    AKTjetMass1 -> Fill(AKTjetpair1Mass);
                    	    AKTjetMass2 -> Fill(AKTjetpair2Mass);
                            //GenAKTMass2->Fill(GenJetMass);
                            AKTGenMass1Comp -> Fill(AKTGenMass1diff);
                            AKTGenMass2Comp -> Fill(AKTGenMass2diff);		
			    AKT1jet1.SetPtEtaPhiM(AKTpt[0], AKTeta[0], AKTphi[0], AKTmass[0]);
			    AKT1jet2.SetPtEtaPhiM(AKTpt[1], AKTeta[1], AKTphi[1], AKTmass[1]);
			    AKT2jet1.SetPtEtaPhiM(AKTpt[2], AKTeta[2], AKTphi[2], AKTmass[2]);
			    AKT2jet2.SetPtEtaPhiM(AKTpt[3], AKTeta[3], AKTphi[3], AKTmass[3]);

			    AKTh1 = AKT1jet1 + AKT1jet2;
			    AKTh2 = AKT2jet1 + AKT2jet2;
			    TLorentzVector AKTdiH = AKTh1 + AKTh2;
                
			    //sort pt
            
			    for (int i = 0; i < 3; i++) {
                   		for (int j = i + 1; j < 4; j++) {
                    		    if (AKTpt[j] > AKTpt[i]) {
                        	        swap(AKTpt[j], AKTpt[i]);
                        		swap(Genpt[j], Genpt[i]);
                        		swap(deltaR[j], deltaR[i]);
                        		swap(AKTeta[j], AKTeta[i]);
					swap(AKTphi[j], AKTphi[i]);
					swap(AKTmass[j], AKTmass[i]);
					swap(AKTBTag[j], AKTBTag[i]);
				    }
				}
			    }
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
			    BDThiggs1pt = AKTh1.Pt();
			    BDThiggs1eta = AKTh1.Eta();
			    BDThiggs1phi = AKTh1.Phi();
			    BDThiggs1invm = AKTh1.Mag();
			    //BDThiggs1invm = AKTjetpair1Mass;

			    BDThiggs2pt = AKTh2.Pt();
			    BDThiggs2eta = AKTh2.Eta();
			    BDThiggs2phi = AKTh2.Phi();
			    BDThiggs2invm = AKTh2.Mag();
			    //BDThiggs2invm = AKTjetpair2Mass;


			    BDThiggsdeltaEta = TMath::Abs(BDThiggs1eta - BDThiggs2eta);
			    BDThiggsdeltaPhi = TMath::Abs(BDThiggs1phi - BDThiggs2phi);

			    BDTdihiggspt = AKTdiH.Pt();
			    BDTdihiggseta = AKTdiH.Eta();
			    BDTdihiggsphi = AKTdiH.Phi();
			    BDTdihiggsinvm = AKTdiH.Mag();
/*
			    if (AKTBTag[0] >=4){
				BDTjet1BTag1 = 1;
			    } else {
				BDTjet1BTag1 = 0;
			    }

			    if (AKTBTag[1] >=4){
				BDTjet1BTag2 = 1;
			    } else {
				BDTjet1BTag2 = 0;
			    }

			    if (AKTBTag[2] >=4){
				BDTjet2BTag1 = 1;
			    } else {
				BDTjet2BTag1 = 0;
			    }

			    if (AKTBTag[3] >=4){
				BDTjet2BTag2 = 1;
			    } else {
				BDTjet2BTag2 = 0;
			    }

			    NBTag = BDTjet1BTag1 + BDTjet1BTag2 + BDTjet2BTag1 + BDTjet2BTag2;
*/
			    DiHiggsMETDeltaTheta = TMath::Abs(2*TMath::ATan(TMath::Exp(-BDTdihiggseta))-2*TMath::ATan(TMath::Exp(-BDTMissingETEta)));

			    tree_BDT -> Fill();


			//}
	
                    //}
                //}
            //}
            }
            //AKTjetMass1->Fill(AKTjetpair1Mass);
            //AKTjetMass2->Fill(AKTjetpair2Mass);

        }

        //GenUncutMass2->Fill(GenJetMass);
    }
    
    //Fitting and plotting
    /*
    cout << endl << "Results of pairing algorithm:" << endl;
    TCanvas *mycanvas = new TCanvas("mycanvas", "My Canvas", 800, 600);
    TF1 *jetpair1fit = new TF1("jetpair1fit", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)", 25, 600);
    TF1 *jetpair2fit = new TF1("jetpair2fit", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)+expo(5)", 25, 600);
    TF1 *fSignal = new TF1("fSignal", "gaus+gaus(3)", 20, 600);
    TF1 *fBackground = new TF1("fBackground", "expo", 20, 600);
    TF1 *f1peak = new TF1("f1peak", "gaus", 20, 600);
    TF1 *f2peak = new TF1("f2peak", "gaus", 20, 600);
    Double_t param[8];

    jetpair2fit -> SetParameters(120, 120, 10, 20, 10, 2, -0.0001);
    jetpair2fit -> SetParLimits(0, 30, 120);
    jetpair2fit -> SetParLimits(1, 100, 140);
    jetpair2fit -> SetParLimits(2, 10, 30);
    jetpair2fit -> SetParLimits(5, 0, 20);
    jetpair2fit -> SetParLimits(6, -1, -0.0001);
    //jetpair2fit->SetParLimits(4,50,109);
    jetpair2fit -> SetParLimits(3, 10, 50);
    jetpair2fit -> SetParLimits(4, 0, 30);

    jetpair1fit -> SetParameters(300, 120, 10, 40, 20);
    jetpair1fit -> SetParLimits(0, 50, 350);
    jetpair1fit -> SetParLimits(1, 100, 140);
    jetpair1fit -> SetParLimits(2, 5, 30);
    jetpair1fit -> SetParLimits(3, 0, 100);
    jetpair1fit -> SetParLimits(4, 5, 40);

    
    cout << endl << "Run gaussian fit for the leading anti-KT jet pair..." << endl;
    AKTjetMass1 -> Fit("jetpair1fit", "R");
    jetpair1fit -> GetParameters(param);
    fSignal -> SetParameters( &param[0]);
    fBackground -> SetParameters( &param[6]);
    f1peak -> SetParameters( &param[0]);
    f2peak -> SetParameters(param[3], param[1], param[4]);
    TH1D *AKTjetMass1Signal = new TH1D( *AKTjetMass1);
    AKTjetMass1Signal -> Sumw2();
    AKTjetMass1Signal -> Add(fBackground, -1);
    AKTjetMass1 -> GetXaxis() -> SetTitle("M [GeV]");
    AKTjetMass1 -> GetYaxis() -> SetTitle("Events");
    AKTjetMass1 -> Draw();
    fBackground -> SetLineColor(4);
    fBackground -> Draw("SAME");
    f1peak -> SetLineColor(3);
    f1peak -> Draw("SAME");
    f2peak -> SetLineColor(6);
    f2peak -> Draw("SAME");
    mycanvas -> SaveAs("AKTjetpair1Mass.png");

    cout << endl << "Run gaussian fit for the sub-leading anti-KT jet pair..." << endl;
    AKTjetMass2 -> Fit("jetpair2fit", "R");
    jetpair2fit -> GetParameters(param);
    fSignal -> SetParameters( &param[0]);
    fBackground -> SetParameters( &param[5]);
    f1peak -> SetParameters( &param[0]);
    f2peak -> SetParameters(param[3], param[1], param[4]);
    TH1D *AKTjetMass2Signal = new TH1D( *AKTjetMass2);
    AKTjetMass2Signal -> Sumw2();
    AKTjetMass2Signal -> Add(fBackground, -1);
    AKTjetMass2 -> GetXaxis() -> SetTitle("M [GeV]");
    AKTjetMass2 -> GetYaxis() -> SetTitle("Events");
    AKTjetMass2 -> Draw();
    fBackground -> SetLineColor(4);
    fBackground -> Draw("SAME");
    f1peak -> SetLineColor(3);
    f1peak -> Draw("SAME");
    f2peak -> SetLineColor(6);
    f2peak -> Draw("SAME");
    //AKTjetMass2Signal->Draw("SAME"); fSignal->Draw("SAME");
    mycanvas -> SaveAs("AKTjetpair2Mass.png");
  
    output -> cd();
    
    cout << endl << "Saving..." << endl;
    AKTjetMass1 -> Write();
    AKTjetMass2 -> Write();
    AKTGenMass1Comp -> Write();
    AKTGenMass2Comp -> Write();
    AKTGenPt1Comp -> Write();
    AKTGenPt2Comp -> Write();
    AKTGenPt3Comp -> Write();
    AKTGenPt4Comp -> Write();

    AKTjetPT_Theta -> Write();
    GenJetPT_Theta -> Write();

    jet1Reso_Pt -> Write();
    jet2Reso_Pt -> Write();
    jet3Reso_Pt -> Write();
    jet4Reso_Pt -> Write();

    jet1Reso_DeltaR -> Write();
    jet2Reso_DeltaR -> Write();
    jet3Reso_DeltaR -> Write();
    jet4Reso_DeltaR -> Write();

    badjet1DeltaR -> Write();
    badjet2DeltaR -> Write();
    badjet3DeltaR -> Write();
    badjet4DeltaR -> Write();
    alljet1DeltaR -> Write();
    alljet2DeltaR -> Write();
    alljet3DeltaR -> Write();
    alljet4DeltaR -> Write();

    badjet1theta -> Write();
    badjet2theta -> Write();
    badjet3theta -> Write();
    badjet4theta -> Write();

    goodjet1theta -> Write();
    goodjet2theta -> Write();
    goodjet3theta -> Write();
    goodjet4theta -> Write();

    alljet1theta -> Write();
    alljet2theta -> Write();
    alljet3theta -> Write();
    alljet4theta -> Write();

    tree_output -> Write();
    output -> Close();
    */
    file_sig -> Close();
    file_BDT -> cd();
    tree_BDT -> Write();
    //printf("Time taken: %.2fs\n", (Double_t)(clock() - tStart)/CLOCKS_PER_SEC);
    cout << "Create dataset for Boosting Decision Tress training...." << endl;
    cout << "Done" << endl;
}
void trainBDT(const char *BDToutputFileName) {
    TMVA::Tools::Instance();
    TFile* BDToutputFile = TFile::Open(BDToutputFileName, "RECREATE");
    TMVA::Factory *factory = new TMVA::Factory("TMVAClassification",BDToutputFile,"!V:!Silent:Color:Transformations=I:DrawProgressBar:AnalysisType=Classification"); 
    TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");
    
    Double_t signalWeight = 1.0;
    Double_t backgroundWeight=1.0;
         
    TFile* dataFile = new TFile("~/Delphes/delphes_dhiggs_sig+bkg_pairmass.root");
    TTree* sigTree = (TTree*)(dataFile -> Get("tree_BDT_sig"));
    TTree* bkg1Tree = (TTree*)(dataFile -> Get("tree_BDT_bkg1"));
    TTree* bkg2Tree = (TTree*)(dataFile -> Get("tree_BDT_bkg2"));
    TTree* bkg3Tree = (TTree*)(dataFile -> Get("tree_BDT_bkg3"));
	      
    dataloader -> TMVA::DataLoader::AddSignalTree(sigTree,signalWeight);
    //dataloader -> TMVA::DataLoader::AddBackgroundTree(bkg1Tree,backgroundWeight);
    //dataloader -> TMVA::DataLoader::AddBackgroundTree(bkg2Tree,backgroundWeight);
    dataloader -> TMVA::DataLoader::AddBackgroundTree(bkg3Tree,backgroundWeight);

    dataloader -> AddVariable("BDTNjets", 'I');
     
    dataloader -> AddVariable("BDTjet1pt1",'F');
    dataloader -> AddVariable("BDTjet1pt2",'F');
    dataloader -> AddVariable("BDTjet2pt1",'F');
    dataloader -> AddVariable("BDTjet2pt2",'F');
    
    dataloader -> AddVariable("BDTjet1eta1",'F');
    dataloader -> AddVariable("BDTjet1eta2",'F');
    dataloader -> AddVariable("BDTjet2eta1",'F');
    dataloader -> AddVariable("BDTjet2eta2",'F');
    
    dataloader -> AddVariable("BDTjet1phi1",'F');
    dataloader -> AddVariable("BDTjet1phi2",'F');
    dataloader -> AddVariable("BDTjet2phi1",'F');
    dataloader -> AddVariable("BDTjet2phi2",'F');
/*   
    dataloader -> AddVariable("BDTjet1BTag1",'I');
    dataloader -> AddVariable("BDTjet1BTag2",'I');
    dataloader -> AddVariable("BDTjet2BTag1",'I');
    dataloader -> AddVariable("BDTjet2BTag2",'I');
*/
    dataloader -> AddVariable("BDTMissingETMET", 'F');
    dataloader -> AddVariable("BDTMissingETEta", 'F');
    dataloader -> AddVariable("BDTMissingETPhi", 'F');

    //dataloader -> AddVariable("NBTag",'I');

    dataloader -> AddVariable("DiHiggsMETDeltaTheta", 'F');
    
    dataloader -> AddVariable("BDThiggsdeltaEta", 'F');
    dataloader -> AddVariable("BDThiggsdeltaPhi", 'F');

    dataloader -> AddVariable("BDThiggs1pt",'F');
    dataloader -> AddVariable("BDThiggs1eta",'F');
    //dataloader -> AddVariable("BDThiggs1phi",'F');
    dataloader -> AddVariable("BDThiggs1invm",'F');

    dataloader -> AddVariable("BDThiggs2pt",'F');
    dataloader -> AddVariable("BDThiggs2eta",'F');
    //dataloader -> AddVariable("BDThiggs2phi",'F');
    dataloader -> AddVariable("BDThiggs2invm",'F');

    dataloader -> AddVariable("BDTdihiggspt",'F');
    dataloader -> AddVariable("BDTdihiggseta",'F');
    //dataloader -> AddVariable("BDTdihiggsphi",'F');
    dataloader -> AddVariable("BDTdihiggsinvm",'F');

    TCut mycuts = TCut("");
    TCut mycutb = TCut("");

    dataloader -> PrepareTrainingAndTestTree(mycuts,mycutb,"random"); 
    factory -> BookMethod(dataloader, TMVA::Types::kBDT,"BDT_RealAdaBoost","!H:!V:NTrees=400:MinNodeSize=3.5%:BoostType=RealAdaBoost:!UseBaggedBoost:BaggedSampleFraction=0.6:nCuts=20:MaxDepth=3:!DoBoostMonitor");
    //factory -> BookMethod(dataloader, TMVA::Types::kMLP, "MLP_ANN", "!H:!V:NeuronType=sigmoid:VarTransform=N:NCycles=2000:HiddenLayers=16,4:TestRate=5:LearningRate=0.05:UseRegulator" );
    //factory -> BookMethod(dataloader, TMVA::Types::kDL, "DNN", "!V:Architecture=CPU:Layout=SIGMOID|1024,SIGMOID|4:TrainingStrategy=LearningRate=2e-4,MaxEpoch=2000");
    
    factory -> TrainAllMethods();
    factory -> TestAllMethods();
    factory -> EvaluateAllMethods();
    
    //cout << endl << "AUC for MLP_ANN method is: " << factory -> GetROCIntegral(dataloader, "MLP_ANN") << endl;
    //cout << endl << "AUC for DNN method is: " << factory -> GetROCIntegral(dataloader, "DNN") << endl;
    //cout << endl << "AUC for BDT_RealAdaBoost method is: " << factory -> GetROCIntegral(dataloader, "BDT_RealAdaBoost") << endl;
 
    BDToutputFile  ->  Close();
    dataFile  ->  Close();
    //TMVA::TMVAGui(BDToutputFileName);
    //gApplication  ->  Run();
}
void applyBDTindividual(TTree* theTree, 
    TH1F* histBdt,
        TH1F* jetpair1,
	    TH1F* jetpair2){
    TMVA::Tools::Instance();
    gROOT -> ProcessLine("gErrorIgnoreLevel = kWarning;");
    gErrorIgnoreLevel = kWarning;
    TMVA::Reader *reader = new TMVA::Reader("Color:!Silent");

    Float_t BDTNjets;

    Float_t BDTjet1pt1;
    Float_t BDTjet1pt2;
    Float_t BDTjet2pt1;
    Float_t BDTjet2pt2;

    Float_t BDTjet1eta1;
    Float_t BDTjet1eta2;
    Float_t BDTjet2eta1;
    Float_t BDTjet2eta2;

    Float_t BDTjet1phi1;
    Float_t BDTjet1phi2;
    Float_t BDTjet2phi1;
    Float_t BDTjet2phi2;
/*
    Float_t BDTjet1BTag1;
    Float_t BDTjet1BTag2;
    Float_t BDTjet2BTag1;
    Float_t BDTjet2BTag2;
*/
    Float_t BDTMissingETMET;
    Float_t BDTMissingETEta;
    Float_t BDTMissingETPhi;

    Float_t BDThiggsdeltaEta;
    Float_t BDThiggsdeltaPhi;

    Float_t BDThiggs1pt;
    Float_t BDThiggs1eta;
    Float_t BDThiggs1invm;

    Float_t BDThiggs2pt;
    Float_t BDThiggs2eta;
    Float_t BDThiggs2invm;

    Float_t BDTdihiggspt;
    Float_t BDTdihiggseta;
    Float_t BDTdihiggsinvm;

    //Float_t NBTag;
    Float_t DiHiggsMETDeltaTheta;

    Int_t BDTNjets_d;

    Double_t BDTjet1pt1_d;
    Double_t BDTjet1pt2_d;
    Double_t BDTjet2pt1_d;
    Double_t BDTjet2pt2_d;

    Double_t BDTjet1eta1_d;
    Double_t BDTjet1eta2_d;
    Double_t BDTjet2eta1_d;
    Double_t BDTjet2eta2_d;

    Double_t BDTjet1phi1_d;
    Double_t BDTjet1phi2_d;
    Double_t BDTjet2phi1_d;
    Double_t BDTjet2phi2_d;
/*
    Int_t BDTjet1BTag1_d;
    Int_t BDTjet1BTag2_d;
    Int_t BDTjet2BTag1_d;
    Int_t BDTjet2BTag2_d;
*/
    Double_t BDTMissingETMET_d;
    Double_t BDTMissingETEta_d;
    Double_t BDTMissingETPhi_d;

    Double_t BDThiggsdeltaEta_d;
    Double_t BDThiggsdeltaPhi_d;

    Double_t BDThiggs1pt_d;
    Double_t BDThiggs1eta_d;
    Double_t BDThiggs1invm_d;

    Double_t BDThiggs2pt_d;
    Double_t BDThiggs2eta_d;
    Double_t BDThiggs2invm_d;

    Double_t BDTdihiggspt_d;
    Double_t BDTdihiggseta_d;
    Double_t BDTdihiggsinvm_d;

    Int_t NBTag_d;
    Double_t DiHiggsMETDeltaTheta_d;


    reader -> AddVariable("BDTNjets", &BDTNjets);

    reader -> AddVariable("BDTjet1pt1", &BDTjet1pt1);
    reader -> AddVariable("BDTjet1pt2", &BDTjet1pt2);
    reader -> AddVariable("BDTjet2pt1", &BDTjet2pt1);
    reader -> AddVariable("BDTjet2pt2", &BDTjet2pt2);

    reader -> AddVariable("BDTjet1eta1", &BDTjet1eta1);
    reader -> AddVariable("BDTjet1eta2", &BDTjet1eta2);
    reader -> AddVariable("BDTjet2eta1", &BDTjet2eta1);
    reader -> AddVariable("BDTjet2eta2", &BDTjet2eta2);

    reader -> AddVariable("BDTjet1phi1", &BDTjet1phi1);
    reader -> AddVariable("BDTjet1phi2", &BDTjet1phi2);
    reader -> AddVariable("BDTjet2phi1", &BDTjet2phi1);
    reader -> AddVariable("BDTjet2phi2", &BDTjet2phi2);
/*
    reader -> AddVariable("BDTjet1BTag1", &BDTjet1BTag1);
    reader -> AddVariable("BDTjet1BTag2", &BDTjet1BTag2);
    reader -> AddVariable("BDTjet2BTag1", &BDTjet2BTag1);
    reader -> AddVariable("BDTjet2BTag2", &BDTjet2BTag2);
*/
    reader -> AddVariable("BDTMissingETMET", &BDTMissingETMET);
    reader -> AddVariable("BDTMissingETEta", &BDTMissingETEta);
    reader -> AddVariable("BDTMissingETPhi", &BDTMissingETPhi);

    //reader -> AddVariable("NBTag", &NBTag);

    reader -> AddVariable("DiHiggsMETDeltaTheta", &DiHiggsMETDeltaTheta);

    reader -> AddVariable("BDThiggsdeltaEta", &BDThiggsdeltaEta);
    reader -> AddVariable("BDThiggsdeltaPhi", &BDThiggsdeltaPhi);

    reader -> AddVariable("BDThiggs1pt", &BDThiggs1pt);
    reader -> AddVariable("BDThiggs1eta", &BDThiggs1eta);
    reader -> AddVariable("BDThiggs1invm", &BDThiggs1invm);

    reader -> AddVariable("BDThiggs2pt", &BDThiggs2pt);
    reader -> AddVariable("BDThiggs2eta", &BDThiggs2eta);
    reader -> AddVariable("BDThiggs2invm", &BDThiggs2invm);

    reader -> AddVariable("BDTdihiggspt", &BDTdihiggspt);
    reader -> AddVariable("BDTdihiggseta", &BDTdihiggseta);
    reader -> AddVariable("BDTdihiggsinvm", &BDTdihiggsinvm);

    reader -> BookMVA("BDT_RealAdaBoost","dataset/weights/TMVAClassification_BDT_RealAdaBoost.weights.xml");
      
    TLeaf *Njets = theTree -> GetLeaf("BDTNjets");

    TLeaf *jet1pt1 = theTree -> GetLeaf("BDTjet1pt1");
    TLeaf *jet1pt2 = theTree -> GetLeaf("BDTjet1pt2");
    TLeaf *jet2pt1 = theTree -> GetLeaf("BDTjet2pt1");
    TLeaf *jet2pt2 = theTree -> GetLeaf("BDTjet2pt2");

    TLeaf *jet1eta1 = theTree -> GetLeaf("BDTjet1eta1");
    TLeaf *jet1eta2 = theTree -> GetLeaf("BDTjet1eta2");
    TLeaf *jet2eta1 = theTree -> GetLeaf("BDTjet2eta1");
    TLeaf *jet2eta2 = theTree -> GetLeaf("BDTjet2eta2");

    TLeaf *jet1phi1 = theTree -> GetLeaf("BDTjet1phi1");
    TLeaf *jet1phi2 = theTree -> GetLeaf("BDTjet1phi2");
    TLeaf *jet2phi1 = theTree -> GetLeaf("BDTjet2phi1");
    TLeaf *jet2phi2 = theTree -> GetLeaf("BDTjet2phi2");
/*
    TLeaf *jet1BTag1 = theTree -> GetLeaf("BDTjet1BTag1");
    TLeaf *jet1BTag2 = theTree -> GetLeaf("BDTjet1BTag2");
    TLeaf *jet2BTag1 = theTree -> GetLeaf("BDTjet2BTag1");
    TLeaf *jet2BTag2 = theTree -> GetLeaf("BDTjet2BTag2");
*/
    //TLeaf *NBTagLeaf = theTree -> GetLeaf("NBTag");

    TLeaf *MissingETMET = theTree -> GetLeaf("BDTMissingETMET");
    TLeaf *MissingETEta = theTree -> GetLeaf("BDTMissingETEta");
    TLeaf *MissingETPhi = theTree -> GetLeaf("BDTMissingETPhi");

    TLeaf *HiggsMetDeltaTheta = theTree -> GetLeaf("DiHiggsMETDeltaTheta");

    TLeaf *higgsdeltaEta = theTree -> GetLeaf("BDThiggsdeltaEta");
    TLeaf *higgsdeltaPhi = theTree -> GetLeaf("BDThiggsdeltaPhi");

    TLeaf *higgs1pt = theTree -> GetLeaf("BDThiggs1pt");
    TLeaf *higgs1eta = theTree -> GetLeaf("BDThiggs1eta");
    TLeaf *higgs1invm = theTree -> GetLeaf("BDThiggs1invm");

    TLeaf *higgs2pt = theTree -> GetLeaf("BDThiggs2pt");
    TLeaf *higgs2eta = theTree -> GetLeaf("BDThiggs2eta");
    TLeaf *higgs2invm = theTree -> GetLeaf("BDThiggs2invm");

    TLeaf *dihiggspt = theTree -> GetLeaf("BDTdihiggspt");
    TLeaf *dihiggseta = theTree -> GetLeaf("BDTdihiggseta");
    TLeaf *dihiggsinvm = theTree -> GetLeaf("BDTdihiggsinvm");

    for (Long64_t entry=0; entry < theTree->GetEntries(); entry++) {
        if (entry%1000 == 0) std::cout << "--- ... Processing event: " << entry << std::endl;
        theTree -> GetEntry(entry);

        BDTNjets = Njets -> GetValue();

        BDTjet1pt1 =  jet1pt1 -> GetValue();
        BDTjet1pt2 =  jet1pt2 -> GetValue();
        BDTjet2pt1 =  jet2pt1 -> GetValue();
        BDTjet2pt2 =  jet2pt2 -> GetValue();

        BDTjet1eta1 =  jet1eta1 -> GetValue();
        BDTjet1eta2 =  jet1eta2 -> GetValue();
        BDTjet2eta1 =  jet2eta1 -> GetValue();
        BDTjet2eta2 =  jet2eta2 -> GetValue();

        BDTjet1phi1 =  jet1phi1 -> GetValue();
        BDTjet1phi2 =  jet1phi2 -> GetValue();
        BDTjet2phi1 =  jet2phi1 -> GetValue();
        BDTjet2phi2 =  jet2phi2 -> GetValue();
/*
        BDTjet1BTag1 = jet1BTag1 -> GetValue();
        BDTjet1BTag2 = jet1BTag2 -> GetValue();
        BDTjet2BTag1 = jet2BTag1 -> GetValue();
        BDTjet2BTag2 = jet2BTag2 -> GetValue();
*/
        BDTMissingETMET = MissingETMET -> GetValue();
        BDTMissingETEta = MissingETEta -> GetValue();
        BDTMissingETPhi = MissingETPhi -> GetValue();

        BDThiggsdeltaEta = higgsdeltaEta -> GetValue();
        BDThiggsdeltaPhi = higgsdeltaPhi -> GetValue();

        BDThiggs1pt =  higgs1pt -> GetValue();
        BDThiggs1eta = higgs1eta -> GetValue();
        BDThiggs1invm = higgs1invm -> GetValue();

        BDThiggs2pt = higgs2pt -> GetValue();
        BDThiggs2eta = higgs2eta -> GetValue();
        BDThiggs2invm = higgs2invm -> GetValue();

        BDTdihiggspt = dihiggspt -> GetValue();
        BDTdihiggseta = dihiggseta -> GetValue();
        BDTdihiggsinvm = dihiggsinvm -> GetValue();

        //NBTag =  NBTagLeaf -> GetValue();
        DiHiggsMETDeltaTheta = HiggsMetDeltaTheta -> GetValue();

        histBdt -> Fill(reader -> EvaluateMVA("BDT_RealAdaBoost"));
        if (reader->EvaluateMVA("BDT_RealAdaBoost") > 0) {
            jetpair1 -> Fill(BDThiggs1invm);
            jetpair2 -> Fill(BDThiggs2invm);
        }
    }
}

void applyBDT(const char *BDTApplyOutputFileName) {
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    TFile* dataFile = new TFile("~/Delphes/delphes_dhiggs_sig+bkg_pairmass.root");
    TTree* sigTree = (TTree*) (dataFile -> Get("tree_BDT_sig"));
    TTree* bkg1Tree = (TTree*) (dataFile -> Get("tree_BDT_bkg1"));
    TTree* bkg2Tree = (TTree*) (dataFile -> Get("tree_BDT_bkg2"));
    TTree* bkg3Tree = (TTree*) (dataFile -> Get("tree_BDT_bkg3"));
    UInt_t nbin = 50;
    TH1F *histBdt_sig = new TH1F("MVA_BDT_sig", "MVA_BDT_sig", nbin, 0.35, 0.65);
    TH1F *histBdt_bkg1 = new TH1F("MVA_BDT_bkg1", "MVA_BDT_bkg1", nbin, 0.35, 0.65);
    TH1F *histBdt_bkg2 = new TH1F("MVA_BDT_bkg2", "MVA_BDT_bkg2", nbin, 0.35, 0.65);
    TH1F *histBdt_bkg3 = new TH1F("MVA_BDT_bkg3", "MVA_BDT_bkg3", nbin, 0.35, 0.65);

    Int_t SigEntries = sigTree -> GetEntries();
    Int_t bkg1Entries = bkg1Tree -> GetEntries();
    Int_t bkg2Entries = bkg2Tree -> GetEntries();
    Int_t bkg3Entries = bkg3Tree -> GetEntries();

    THStack* JetPair1BDT = new THStack("JetPair1BDT", "");
    THStack* JetPair2BDT = new THStack("JetPair2BDT", "");
    TH1F *AKTjetMass1_sigBDT = new TH1F("AKTjetMass1_sigBDT", "Anti_KTjet leading jets pair invariant mass", 50, 80, 160);
    TH1F *AKTjetMass2_sigBDT = new TH1F("AKTjetMass2_sigBDT", "Anti_KTjet sub-leading jets pair invariant mass", 50, 80, 160);

    TH1F *AKTjetMass1_bkg1BDT = new TH1F("AKTjetMass1_bkg1BDT", "Anti_KTjet leading jets pair invariant mass", 50, 80, 160);
    TH1F *AKTjetMass2_bkg1BDT = new TH1F("AKTjetMass2_bkg1BDT", "Anti_KTjet sub-leading jets pair invariant mass", 50, 80, 160);

    TH1F *AKTjetMass1_bkg2BDT = new TH1F("AKTjetMass1_bkg2BDT", "Anti_KTjet leading jets pair invariant mass", 50, 80, 160);
    TH1F *AKTjetMass2_bkg2BDT = new TH1F("AKTjetMass2_bkg2BDT", "Anti_KTjet sub-leading jets pair invariant mass", 50, 80, 160);

    TH1F *AKTjetMass1_bkg3BDT = new TH1F("AKTjetMass1_bkg3BDT", "Anti_KTjet leading jets pair invariant mass", 50, 80, 160);
    TH1F *AKTjetMass2_bkg3BDT = new TH1F("AKTjetMass2_bkg3BDT", "Anti_KTjet sub-leading jets pair invariant mass", 50, 80, 160);

    Double_t weight1 = 100*0.0008201/(0.0008201 + 0.0009227 + 0.003743 + 0.03168);
    Double_t weight2 = 100*0.0009227/(0.0008201 + 0.0009227 + 0.003743 + 0.03168);
    Double_t weight3 = 100*0.003743/(0.0008201 + 0.0009227 + 0.003743 + 0.03168);
    Double_t weight4 = 100*0.03168/(0.0008201 + 0.0009227 + 0.003743 + 0.03168);

    Double_t SigStrength = SigEntries * weight1/TMath::Sqrt(SigEntries * weight1 + bkg1Entries * weight2 + bkg2Entries * weight3 + bkg3Entries * weight4);

    cout << SigStrength << endl;


    applyBDTindividual(sigTree, histBdt_sig, AKTjetMass1_sigBDT, AKTjetMass2_sigBDT);
    applyBDTindividual(bkg1Tree, histBdt_bkg1, AKTjetMass1_bkg1BDT, AKTjetMass2_bkg1BDT);
    applyBDTindividual(bkg2Tree, histBdt_bkg2, AKTjetMass1_bkg2BDT, AKTjetMass2_bkg2BDT);
    applyBDTindividual(bkg3Tree, histBdt_bkg3, AKTjetMass1_bkg3BDT, AKTjetMass2_bkg3BDT);

    //Stack histo for invm
    AKTjetMass1_sigBDT -> Scale(weight1);
    //AKTjetMass1_sigBDT -> SetFillColor(kPink-1);
    AKTjetMass1_sigBDT -> SetLineColor(kPink-1);
    AKTjetMass1_sigBDT -> SetLineWidth(4);
    //JetPair1BDT -> Add(AKTjetMass1_sigBDT);
    AKTjetMass1_bkg1BDT -> Scale(weight2);
    AKTjetMass1_bkg1BDT-> SetFillColor(kAzure-1);
    JetPair1BDT -> Add(AKTjetMass1_bkg1BDT);
    AKTjetMass1_bkg2BDT -> Scale(weight3);
    AKTjetMass1_bkg2BDT-> SetFillColor(kAzure+1);
    JetPair1BDT -> Add(AKTjetMass1_bkg2BDT);
    AKTjetMass1_bkg3BDT -> Scale(weight4);
    AKTjetMass1_bkg3BDT-> SetFillColor(kAzure+8);
    JetPair1BDT -> Add(AKTjetMass1_bkg3BDT);
    AKTjetMass2_sigBDT -> Scale(weight1);
    //AKTjetMass2_sigBDT -> SetFillColor(kPink-1);
    AKTjetMass2_sigBDT -> SetLineColor(kPink-1);
    AKTjetMass2_sigBDT -> SetLineWidth(4);
    //JetPair2BDT -> Add(AKTjetMass2_sigBDT);
    AKTjetMass2_bkg1BDT -> Scale(weight2);
    AKTjetMass2_bkg1BDT-> SetFillColor(kAzure-1);
    JetPair2BDT -> Add(AKTjetMass2_bkg1BDT);
    AKTjetMass2_bkg2BDT -> Scale(weight3);
    AKTjetMass2_bkg2BDT-> SetFillColor(kAzure+1);
    JetPair2BDT -> Add(AKTjetMass2_bkg2BDT);
    AKTjetMass2_bkg3BDT -> Scale(weight4);
    AKTjetMass2_bkg3BDT-> SetFillColor(kAzure+8);
    JetPair2BDT -> Add(AKTjetMass2_bkg3BDT);

    TF1 *fSignal = new TF1("fSignal", "[0]*exp(-0.5*((x-[1])/[2])^2)", 80, 160);
    TF1 *f1Signal = new TF1("Signal1", "[0]*exp(-0.5*((x-[1])/[2])^2)", 80, 160);
    TF1 *f2Bkg = new TF1("Bkg2", "[0]*exp(-0.5*((x-[1])/[2])^2)", 80, 160);
    TF1 *f3Bkg = new TF1("Bkg3", "[0]*exp(-0.5*((x-[1])/[2])^2)", 80, 160);
    //TF1 *fAll = new TF1("All", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)+[6]*exp(-0.5*((x-[7])/[8])^2)", 80, 160);

    Double_t param1[3];
    Double_t param2[3];
    Double_t param3[3];

    fSignal -> SetParameters(70, 120, 40);

    TCanvas *totalcanvas = new TCanvas("totalcanvas", "Canvas", 300, 300, 1200, 1200);
    gStyle->SetHistMinimumZero();
    totalcanvas -> SetLogy();
    totalcanvas -> SetWindowSize(1500,1500);
    totalcanvas -> SetCanvasSize(1200,1200);
    JetPair1BDT -> Draw("HIST");
    AKTjetMass1_sigBDT -> Draw("sameHIST");

    AKTjetMass1_sigBDT -> Fit("fSignal", "R");
    fSignal -> GetParameters(param1);
    f1Signal -> SetParameters(&param1[0]);
    f1Signal -> SetLineColor(2);
    f1Signal -> Draw("sameHIST");

/*
    AKTjetMass1_bkg2BDT -> Fit("fSignal", "R");
    fSignal -> GetParameters(param2);

    AKTjetMass1_bkg3BDT -> Fit("fSignal", "R");
    fSignal -> GetParameters(param3);
*/
/*
    fAll -> SetParameters(param1[0],param1[1],param1[2],param2[0],param2[1],param2[2],param3[0],param3[1],param3[2]);
    fAll -> SetLineColor(kAzure+8);
    fAll -> Draw("same");
*/
    TLegend *legend = new TLegend(0.65, 0.65, 0.875, 0.85);
    legend -> AddEntry(AKTjetMass1_sigBDT, "signal HH data", "l");
    legend -> AddEntry(AKTjetMass1_bkg1BDT, "QCD", "f");
    legend -> AddEntry(AKTjetMass1_bkg2BDT, "single Higgs", "f");
    legend -> AddEntry(AKTjetMass1_bkg3BDT, "Z + bb", "f");
    legend -> AddEntry(f1Signal, "Sig HH fit", "L");
    //legend -> AddEntry(fAll, "Sig+Bkg fit", "L");
    legend -> SetBorderSize(0);
    legend -> Draw();
    totalcanvas -> cd();
    JetPair1BDT -> GetXaxis() -> SetTitle("m_{H_1} [GeV]");
    JetPair1BDT -> SetMaximum(1e5);
    JetPair1BDT -> SetMinimum(4);
    JetPair1BDT -> GetYaxis() -> SetTitle("events");
    TLatex title1(95, 2.4e5, "#bf{#font[72]{Muon Collider Simulation (Delphes)}}");
    title1.SetTextSize(0.04);
    TLatex latexUp1(85, 5e4, "\\sqrt{s} = 3 TeV}");
    TLatex latexDown1(88, 2e4, "#font[52]{L} #font[42]{=} #font[52]{1 ab^{-1}}");
    string init("#font[52]{#frac{S}{#sqrt{S+B}}} #font[42]{=} #font[52]{");
    string add = to_string(SigStrength);
    string end("}");
    init = init + add + end;
    const char * latex = init.c_str();
    TLatex latexSigStrength(83, 8e3, latex);
    latexUp1.SetTextSize(0.04);
    latexDown1.SetTextSize(0.04);
    latexSigStrength.SetTextSize(0.04);
    latexUp1.Draw("same");
    latexDown1.Draw("same");
    latexSigStrength.Draw("same");
    title1.Draw("same");
    totalcanvas -> SaveAs("JetPairs1_sig+bkgBDT.png");

    JetPair2BDT -> Draw("HIST");
    AKTjetMass2_sigBDT -> Draw("sameHIST");
    TLegend *legend2 = new TLegend(0.65, 0.65, 0.875, 0.85);
    legend2 -> AddEntry(AKTjetMass2_sigBDT, "signal HH", "l");
    legend2 -> AddEntry(AKTjetMass2_bkg1BDT, "QCD", "f");
    legend2 -> AddEntry(AKTjetMass2_bkg2BDT, "single Higgs ", "f");
    legend2 -> AddEntry(AKTjetMass2_bkg3BDT, "Z + bb", "f");
    legend2 -> SetBorderSize(0);
    legend2 -> Draw();
    JetPair2BDT -> GetXaxis() -> SetTitle("m_{H_2} [GeV]");
    JetPair2BDT -> GetYaxis() -> SetTitle("events");
    JetPair2BDT -> SetMaximum(1e5);
    JetPair2BDT -> SetMinimum(1);
    TLatex title2(95, 2.5e5, "#bf{#font[72]{Muon Collider Simulation (Delphes)}}");
    title2.SetTextSize(0.04);
    TLatex latexUp2(85, 5e4, "\\sqrt{s} = 3 TeV}");
    TLatex latexDown2(90, 1.75e4, "#font[52]{L} #font[42]{=} #font[52]{1 ab^{-1}}");
    latexUp2.SetTextSize(0.04);
    latexDown2.SetTextSize(0.04);
    latexUp2.Draw("same");
    latexDown2.Draw("same");
    title2.Draw("same");
    totalcanvas -> SaveAs("JetPairs2_sig+bkgBDT.png");

//Stack histo for BDT

    THStack* BDTs_b = new THStack("BDTs_b", "");

    histBdt_sig -> Scale(weight1);
    //AKTjetMass1_sigBDT -> SetFillColor(kPink-1);
    histBdt_sig -> SetLineColor(kPink-1);
    histBdt_sig -> SetLineWidth(4);
    //JetPair1BDT -> Add(AKTjetMass1_sigBDT);
    histBdt_bkg1 -> Scale(weight2);
    histBdt_bkg1 -> SetFillColor(kAzure-1);
    BDTs_b -> Add(histBdt_bkg1);
    histBdt_bkg2 -> Scale(weight3);
    histBdt_bkg2 -> SetFillColor(kAzure+1);
    BDTs_b -> Add(histBdt_bkg2);
    histBdt_bkg3  -> Scale(weight4);
    histBdt_bkg3 -> SetFillColor(kAzure+8);
    BDTs_b -> Add(histBdt_bkg3);

    BDTs_b -> Draw("HIST");
    BDTs_b -> SetMaximum(1e4);
    histBdt_sig -> Draw("sameHIST");
    TLegend *legendBdt = new TLegend(0.65, 0.65, 0.875, 0.85);
    legendBdt -> AddEntry(histBdt_sig, "signal HH", "l");
    legendBdt -> AddEntry(histBdt_bkg1, "QCD", "f");
    legendBdt -> AddEntry(histBdt_bkg2, "single Higgs", "f");
    legendBdt -> AddEntry(histBdt_bkg3, "Z + bb", "f");
    legendBdt -> SetBorderSize(0);
    legendBdt -> Draw();
    BDTs_b -> GetXaxis() -> SetTitle("BDT_{Z+bb} score");
    BDTs_b -> GetYaxis() -> SetTitle("events");
    TLatex titleBdt(0.41, 2.15e4, "#bf{#font[72]{Muon Collider Simulation (Delphes)}}");
    titleBdt.SetTextSize(0.04);
    TLatex latexUpBdt(0.365, 5000, "\\sqrt{s} = 3 TeV}");
    TLatex latexDownBdt(0.375, 3000, "#font[52]{L = 1 ab^{-1}}");
    latexUpBdt.SetTextSize(0.04);
    latexDownBdt.SetTextSize(0.04);
    latexUpBdt.Draw("same");
    latexDownBdt.Draw("same");
    titleBdt.Draw("same");
    totalcanvas -> SaveAs("BDT_sig+bkgBDT.png");

    TFile* BDTApplyOutputFile = TFile::Open(BDTApplyOutputFileName, "RECREATE");
    histBdt_sig -> Write();
    histBdt_bkg1 -> Write();
    histBdt_bkg2 -> Write();
    histBdt_bkg3 -> Write();

    BDTApplyOutputFile -> Close();
    dataFile -> Close();
}

void Delphes_HH(const char *inputSigFile,
    const char *inputBkg1File,
    	const char *inputBkg2File,
	    const char *inputBkg3File,
	        const char *outputFile,
		    const char *inputFileForJES,
			const char *outputFileForJES) {
    struct timespec begin, end;
    clock_gettime(CLOCK_REALTIME, &begin);
    //Initiation
    gSystem -> Load("libDelphes.so");
    gROOT -> ProcessLine("gErrorIgnoreLevel = kWarning;");
    gErrorIgnoreLevel = kWarning;

    //Calibration
    Double_t JER[10][10];
    cout << endl << "Calling Calibration algorithm..." << endl;
    Calibration(inputFileForJES, outputFileForJES, JER);

    TFile *output = new TFile(outputFile, "RECREATE");
    TTree *tree_output = new TTree("tree_output", "Delphes");
    TTree *tree_BDT_sig = new TTree("tree_BDT_sig", "sig");
    TTree *tree_BDT_bkg1 = new TTree("tree_BDT_bkg1", "bkg1");
    TTree *tree_BDT_bkg2= new TTree("tree_BDT_bkg2", "bkg2");
    TTree *tree_BDT_bkg3 = new TTree("tree_BDT_bkg3", "bkg3");
    
    THStack* JetPair1 = new THStack("JetPair1", "");// "Leading Anti_{kT} jets pair invariant mass");
    THStack* JetPair2 = new THStack("JetPair2", "");//"Sub-Leading Anti_{kT} jets pair invariant mass");

    TH1D *AKTjetMass1_sig = new TH1D("AKTjetMass1_sig", "Anti_KTjet leading jets pair invariant mass", 50, 0, 250);
    TH1D *AKTjetMass2_sig = new TH1D("AKTjetMass2_sig", "Anti_KTjet sub-leading jets pair invariant mass", 50, 0, 250); 

    TH1D *AKTjetMass1_bkg1 = new TH1D("AKTjetMass1_bkg1", "Anti_KTjet leading jets pair invariant mass", 50, 0, 250);
    TH1D *AKTjetMass2_bkg1 = new TH1D("AKTjetMass2_bkg1", "Anti_KTjet sub-leading jets pair invariant mass", 50, 0, 250); 


    TH1D *AKTjetMass1_bkg2 = new TH1D("AKTjetMass1_bkg2", "Anti_KTjet leading jets pair invariant mass", 50, 0, 250);
    TH1D *AKTjetMass2_bkg2 = new TH1D("AKTjetMass2_bkg2", "Anti_KTjet sub-leading jets pair invariant mass", 50, 0, 250); 


    TH1D *AKTjetMass1_bkg3 = new TH1D("AKTjetMass1_bkg3", "Anti_KTjet leading jets pair invariant mass", 50, 0, 250);
    TH1D *AKTjetMass2_bkg3 = new TH1D("AKTjetMass2_bkg3", "Anti_KTjet sub-leading jets pair invariant mass", 50, 0, 250); 

    Pairing_w_JES(inputSigFile, JER, AKTjetMass1_sig, AKTjetMass2_sig, tree_BDT_sig, output); 
    Pairing_w_JES(inputBkg1File, JER, AKTjetMass1_bkg1, AKTjetMass2_bkg1, tree_BDT_bkg1, output); 
    Pairing_w_JES(inputBkg2File, JER, AKTjetMass1_bkg2, AKTjetMass2_bkg2, tree_BDT_bkg2, output); 
    Pairing_w_JES(inputBkg3File, JER, AKTjetMass1_bkg3, AKTjetMass2_bkg3, tree_BDT_bkg3, output); 

    Int_t N1sig = AKTjetMass1_sig -> GetEntries(); 
    Int_t N1bkg1 = AKTjetMass1_bkg1 -> GetEntries(); 
    Int_t N1bkg2 = AKTjetMass1_bkg2 -> GetEntries(); 
    Int_t N1bkg3 = AKTjetMass1_bkg3 -> GetEntries(); 
    Int_t N2sig = AKTjetMass2_sig -> GetEntries(); 
    Int_t N2bkg1 = AKTjetMass2_bkg1 -> GetEntries(); 
    Int_t N2bkg2 = AKTjetMass2_bkg2 -> GetEntries(); 
    Int_t N2bkg3 = AKTjetMass2_bkg3 -> GetEntries(); 

    Double_t weight1 = 0.0008201/(0.0008201 + 0.0009227 + 0.003743 + 0.03168);
    Double_t weight2 = 0.0009227/(0.0008201 + 0.0009227 + 0.003743 + 0.03168);
    Double_t weight3 = 0.003743/(0.0008201 + 0.0009227 + 0.003743 + 0.03168);
    Double_t weight4 = 0.03168/(0.0008201 + 0.0009227 + 0.003743 + 0.03168);

    AKTjetMass1_sig -> Scale(weight1);
    AKTjetMass1_sig -> SetFillColor(kPink-1);
    JetPair1 -> Add(AKTjetMass1_sig);
    AKTjetMass1_bkg1 -> Scale(weight2);
    AKTjetMass1_bkg1-> SetFillColor(kAzure-1);
    JetPair1 -> Add(AKTjetMass1_bkg1);
    AKTjetMass1_bkg2 -> Scale(weight3);
    AKTjetMass1_bkg2-> SetFillColor(kAzure+1);
    JetPair1 -> Add(AKTjetMass1_bkg2);
    AKTjetMass1_bkg3 -> Scale(weight4);
    AKTjetMass1_bkg3-> SetFillColor(kAzure+8);
    JetPair1 -> Add(AKTjetMass1_bkg3);
 
    AKTjetMass2_sig -> Scale(weight1);
    AKTjetMass2_sig -> SetFillColor(kPink-1);
    JetPair2 -> Add(AKTjetMass2_sig);
    AKTjetMass2_bkg1 -> Scale(weight2);
    AKTjetMass2_bkg1-> SetFillColor(kAzure-1);
    JetPair2 -> Add(AKTjetMass2_bkg1);
    AKTjetMass2_bkg2 -> Scale(weight3);
    AKTjetMass2_bkg2-> SetFillColor(kAzure+1);
    JetPair2 -> Add(AKTjetMass2_bkg2);
    AKTjetMass2_bkg3 -> Scale(weight4);
    AKTjetMass2_bkg3-> SetFillColor(kAzure+8);
    JetPair2 -> Add(AKTjetMass2_bkg3);

    TCanvas *totalcanvas = new TCanvas("totalcanvas", "Canvas", 10, 10, 800, 800);
    //totalcanvas -> SetLogy();
    JetPair1 -> Draw("HIST");
    TLegend *legend = new TLegend(0.6, 0.65, 0.85, 0.85);
    legend -> AddEntry(AKTjetMass1_sig, "  Signal #minus #nu#bar{#nu}HH", "f");
    legend -> AddEntry(AKTjetMass1_bkg1, "  Background #minus #nu#bar{#nu}b#bar{b}b#bar{b}", "f");
    legend -> AddEntry(AKTjetMass1_bkg2, "  Background #minus #nu#bar{#nu}b#bar{b}H", "f");
    legend -> AddEntry(AKTjetMass1_bkg3, "  Background #minus #nu#bar{#nu}b#bar{b}Z", "f");
    legend -> SetBorderSize(0);
    legend -> Draw();
    JetPair1 -> GetXaxis() -> SetTitle("m_{H_1} [GeV]");
    JetPair1 -> GetYaxis() -> SetTitle("Events");
    TLatex title1(30, 315, "#bf{#font[72]{Muon Collider Simulation (Delphes)}}");
    title1.SetTextSize(0.045); 
    TLatex latex1(15, 270, "#it{#sqrt{s}  = 3 TeV}");
    latex1.SetTextSize(0.045); 
    latex1.Draw("same");
    title1.Draw("same");
    totalcanvas -> SaveAs("JetPairs1_sig+bkg.png");

    JetPair2 -> Draw("HIST");
    TLegend *legend2 = new TLegend(0.6, 0.65, 0.85, 0.85);
    legend2 -> AddEntry(AKTjetMass2_sig, "  Signal #minus #nu#bar{#nu}HH", "f");
    legend2 -> AddEntry(AKTjetMass2_bkg1, "  Background #minus #nu#bar{#nu}b#bar{b}b#bar{b}", "f");
    legend2 -> AddEntry(AKTjetMass2_bkg2, "  Background #minus #nu#bar{#nu}b#bar{b}H", "f");
    legend2 -> AddEntry(AKTjetMass2_bkg3, "  Background #minus #nu#bar{#nu}b#bar{b}Z", "f");
    legend2 -> SetBorderSize(0);
    legend2 -> Draw();
    JetPair2 -> GetXaxis() -> SetTitle("m_{H_2} [GeV]");
    JetPair2 -> GetYaxis() -> SetTitle("Events");
    TLatex title2(30, 285, "#bf{#font[72]{Muon Collider Simulation (Delphes)}}");
    title2.SetTextSize(0.045); 
    TLatex latex2(15, 240, "#it{#sqrt{s}  = 3 TeV}");
    latex2.SetTextSize(0.045); 
    latex2.Draw("same");
    title2.Draw("same");
    totalcanvas -> SaveAs("JetPairs2_sig+bkg.png");

    output -> cd();

    JetPair1 -> Write();
    JetPair2 -> Write();

    tree_output -> Write();
    output -> Close();
    //train BDT model
    cout << "Initating TMVA for Multi-Variate Analysis...." << endl;
    trainBDT("delphes_dhiggs_sig+bkg_BDT.root");
    //apply BDT model
    cout << "Start Boosting Decision Trees application..." << endl;
    applyBDT("delphes_dhiggs_sig+bkg_BDTapply.root");

    clock_gettime(CLOCK_REALTIME, &end);
    long seconds = end.tv_sec - begin.tv_sec;
    long nanoseconds = end.tv_nsec - begin.tv_nsec;
    double elapsed = seconds + nanoseconds*1e-9;
    printf("Time taken: %.3f seconds.\n",elapsed);
    cout << "All done! Exit..." << endl;

}
