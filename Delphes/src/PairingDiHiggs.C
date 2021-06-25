//usage: root -l PairingDiHiggs.C\(\"inputfile.root\"\,\"outputfile.root\"\)
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
void JetEnergyFix(Double_t AKTjeteta, Double_t AKTjetPt, Double_t JER[10][10]){
     Double_t AKTjetTheta = 2 * atan(exp(-AKTjeteta));
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
     if (1.8< AKTjetTheta and AKTjetTheta <= 2.1) {
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
     AKTjetPt = AKTjetPt/JER[ThetaGrid][PTGrid];
}

void PairingDiHiggs(const char *inputFile, const char *outputFile){
     gSystem->Load("libDelphes.so");
     //TChain chain("Delphes");
     //chain.Add(inputFile);
     //ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
     TFile *file_sig = new TFile(inputFile);
     TFile *output = new TFile(outputFile, "recreate");
     TTree *tree_sig = (TTree*)file_sig->Get("Delphes");
     TTree *tree_output = new TTree("tree_output","Delphes");

     cout << "Initiating..." <<endl; 
     
     TLeaf *AKTjet_size = tree_sig->GetLeaf("AKTjet_size");
     TLeaf *AKTjet_eta = tree_sig->GetLeaf("AKTjet.Eta");
     TLeaf *AKTjet_phi = tree_sig->GetLeaf("AKTjet.Phi");
     TLeaf *AKTjet_pt = tree_sig->GetLeaf("AKTjet.PT");
     TLeaf *AKTjet_mass = tree_sig->GetLeaf("AKTjet.Mass");
 
     TLeaf *GenJet_size = tree_sig->GetLeaf("GenJet_size");
     TLeaf *GenJet_eta = tree_sig->GetLeaf("GenJet.Eta");
     TLeaf *GenJet_phi = tree_sig->GetLeaf("GenJet.Phi");
     TLeaf *GenJet_pt = tree_sig->GetLeaf("GenJet.PT");
     TLeaf *GenJet_mass = tree_sig->GetLeaf("GenJet.Mass");

     Int_t nEntries = tree_sig->GetEntries();

     TH1D *AKTjetMass1 = new TH1D("AKTjetMass1", "Anti_KTjet leading jets pair invariant mass", 150 , 0, 400); 
     TH1D *AKTjetMass2 = new TH1D("AKTjetMass2", "Anti_KTjet sub-leading jets pair invariant mass", 150 , 0, 400); 
     TH1D *GenAKTMass2 = new TH1D("GenAKTMass2", "GenAKTMass2", 50 , 0, 600); 
     TH1D *AKTGenMass1Comp = new TH1D("AKTGenMass1Comp", "AKTGenMass1Comp", 50 , -1, 2); 
     TH1D *AKTGenMass2Comp = new TH1D("AKTGenMass2Comp", "AKTGenMass2Comp", 50 , -1, 2); 
     TH1D *GenUncutMass2 = new TH1D("GenUncutMass2", "GenUncutMass2", 50 , 0, 600); 
 
     TH1D *AKTGenPt1Comp = new TH1D("AKTGenPt1Comp", "AKTGenPt1Comp", 200 , -1, 5); 
     TH1D *AKTGenPt2Comp = new TH1D("AKTGenPt2Comp", "AKTGenPt2Comp", 200 , -1, 5); 
     TH1D *AKTGenPt3Comp = new TH1D("AKTGenPt3Comp", "AKTGenPt3Comp", 200 , -1, 5); 
     TH1D *AKTGenPt4Comp = new TH1D("AKTGenPt4Comp", "AKTGenPt4Comp", 200 , -1, 5); 

     TH2D *AKTjetPT_Theta = new TH2D("AKTjetPT_Theta", "AKTjetPT_Theta", 30, 0.1482, 3, 30, 0, 400);
     TH2D *GenJetPT_Theta = new TH2D("GenJetPT_Theta", "GenJetPT_Theta", 30, 0.1482, 3, 30, 0, 400);

     TH2D *jet1Reso_Pt = new TH2D("jet1Reso_Pt", "jet1Reso_Pt", 70, 0, 400, 70, -1, 4); 
     TH2D *jet2Reso_Pt = new TH2D("jet2Reso_Pt", "jet2Reso_Pt", 70, 0, 400, 70, -1, 4); 
     TH2D *jet3Reso_Pt = new TH2D("jet3Reso_Pt", "jet3Reso_Pt", 70, 0, 400, 70, -1, 4); 
     TH2D *jet4Reso_Pt = new TH2D("jet4Reso_Pt", "jet4Reso_Pt", 70, 0, 400, 70, -1, 4); 
/*
     TH1D *badjet1eta = new TH1D("badjet1eta", "badjet1eta", 100, -3, 3); 
     TH1D *badjet2eta = new TH1D("badjet2eta", "badjet2eta", 100, -3, 3); 
     TH1D *badjet3eta = new TH1D("badjet3eta", "badjet3eta", 100, -3, 3); 
     TH1D *badjet4eta = new TH1D("badjet4eta", "badjet4eta", 100, -3, 3); 

     TH1D *goodjet1eta = new TH1D("goodjet1eta", "goodjet1eta", 100, -3, 3); 
     TH1D *goodjet2eta = new TH1D("goodjet2eta", "goodjet2eta", 100, -3, 3); 
     TH1D *goodjet3eta = new TH1D("goodjet3eta", "goodjet3eta", 100, -3, 3); 
     TH1D *goodjet4eta = new TH1D("goodjet4eta", "goodjet4eta", 100, -3, 3); 

     TH1D *alljet1eta = new TH1D("alljet1eta", "alljet1eta", 100, -3, 3); 
     TH1D *alljet2eta = new TH1D("alljet2eta", "alljet2eta", 100, -3, 3); 
     TH1D *alljet3eta = new TH1D("alljet3eta", "alljet3eta", 100, -3, 3); 
     TH1D *alljet4eta = new TH1D("alljet4eta", "alljet4eta", 100, -3, 3); 
*/
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
     
     //Calibration histogram
     TH1D *jetPTresponse00 = new TH1D("jetPTreponse00", "jetPTresponse00", 50, 0, 1.3);
     TH1D *jetPTresponse01 = new TH1D("jetPTreponse01", "jetPTresponse01", 50, 0, 1.3);
     TH1D *jetPTresponse02 = new TH1D("jetPTreponse02", "jetPTresponse02", 50, 0, 1.3);
     TH1D *jetPTresponse03 = new TH1D("jetPTreponse03", "jetPTresponse03", 50, 0, 1.3);
     TH1D *jetPTresponse04 = new TH1D("jetPTreponse04", "jetPTresponse04", 50, 0, 1.3);
     TH1D *jetPTresponse05 = new TH1D("jetPTreponse05", "jetPTresponse05", 50, 0, 1.3);
     TH1D *jetPTresponse06 = new TH1D("jetPTreponse06", "jetPTresponse06", 50, 0, 1.3);
     TH1D *jetPTresponse07 = new TH1D("jetPTreponse07", "jetPTresponse07", 50, 0, 1.3);
     TH1D *jetPTresponse08 = new TH1D("jetPTreponse08", "jetPTresponse08", 50, 0, 1.3);
     TH1D *jetPTresponse09 = new TH1D("jetPTreponse09", "jetPTresponse09", 50, 0, 1.3);

     TH1D *jetPTresponse10 = new TH1D("jetPTreponse10", "jetPTresponse10", 50, 0, 1.3);
     TH1D *jetPTresponse11 = new TH1D("jetPTreponse11", "jetPTresponse11", 50, 0, 1.3);
     TH1D *jetPTresponse12 = new TH1D("jetPTreponse12", "jetPTresponse12", 50, 0, 1.3);
     TH1D *jetPTresponse13 = new TH1D("jetPTreponse13", "jetPTresponse13", 50, 0, 1.3);
     TH1D *jetPTresponse14 = new TH1D("jetPTreponse14", "jetPTresponse14", 50, 0, 1.3);
     TH1D *jetPTresponse15 = new TH1D("jetPTreponse15", "jetPTresponse15", 50, 0, 1.3);
     TH1D *jetPTresponse16 = new TH1D("jetPTreponse16", "jetPTresponse16", 50, 0, 1.3);
     TH1D *jetPTresponse17 = new TH1D("jetPTreponse17", "jetPTresponse17", 50, 0, 1.3);
     TH1D *jetPTresponse18 = new TH1D("jetPTreponse18", "jetPTresponse18", 50, 0, 1.3);
     TH1D *jetPTresponse19 = new TH1D("jetPTreponse19", "jetPTresponse19", 50, 0, 1.3);

     TH1D *jetPTresponse20 = new TH1D("jetPTreponse20", "jetPTresponse20", 50, 0, 1.3);
     TH1D *jetPTresponse21 = new TH1D("jetPTreponse21", "jetPTresponse21", 50, 0, 1.3);
     TH1D *jetPTresponse22 = new TH1D("jetPTreponse22", "jetPTresponse22", 50, 0, 1.3);
     TH1D *jetPTresponse23 = new TH1D("jetPTreponse23", "jetPTresponse23", 50, 0, 1.3);
     TH1D *jetPTresponse24 = new TH1D("jetPTreponse24", "jetPTresponse24", 50, 0, 1.3);
     TH1D *jetPTresponse25 = new TH1D("jetPTreponse25", "jetPTresponse25", 50, 0, 1.3);
     TH1D *jetPTresponse26 = new TH1D("jetPTreponse26", "jetPTresponse26", 50, 0, 1.3);
     TH1D *jetPTresponse27 = new TH1D("jetPTreponse27", "jetPTresponse27", 50, 0, 1.3);
     TH1D *jetPTresponse28 = new TH1D("jetPTreponse28", "jetPTresponse28", 50, 0, 1.3);
     TH1D *jetPTresponse29 = new TH1D("jetPTreponse29", "jetPTresponse29", 50, 0, 1.3);


     TH1D *jetPTresponse30 = new TH1D("jetPTreponse30", "jetPTresponse30", 50, 0, 1.3);
     TH1D *jetPTresponse31 = new TH1D("jetPTreponse31", "jetPTresponse31", 50, 0, 1.3);
     TH1D *jetPTresponse32 = new TH1D("jetPTreponse32", "jetPTresponse32", 50, 0, 1.3);
     TH1D *jetPTresponse33 = new TH1D("jetPTreponse33", "jetPTresponse33", 50, 0, 1.3);
     TH1D *jetPTresponse34 = new TH1D("jetPTreponse34", "jetPTresponse34", 50, 0, 1.3);
     TH1D *jetPTresponse35 = new TH1D("jetPTreponse35", "jetPTresponse35", 50, 0, 1.3);
     TH1D *jetPTresponse36 = new TH1D("jetPTreponse36", "jetPTresponse36", 50, 0, 1.3);
     TH1D *jetPTresponse37 = new TH1D("jetPTreponse37", "jetPTresponse37", 50, 0, 1.3);
     TH1D *jetPTresponse38 = new TH1D("jetPTreponse38", "jetPTresponse38", 50, 0, 1.3);
     TH1D *jetPTresponse39 = new TH1D("jetPTreponse39", "jetPTresponse39", 50, 0, 1.3);


     TH1D *jetPTresponse40 = new TH1D("jetPTreponse40", "jetPTresponse40", 50, 0, 1.3);
     TH1D *jetPTresponse41 = new TH1D("jetPTreponse41", "jetPTresponse41", 50, 0, 1.3);
     TH1D *jetPTresponse42 = new TH1D("jetPTreponse42", "jetPTresponse42", 50, 0, 1.3);
     TH1D *jetPTresponse43 = new TH1D("jetPTreponse43", "jetPTresponse43", 50, 0, 1.3);
     TH1D *jetPTresponse44 = new TH1D("jetPTreponse44", "jetPTresponse44", 50, 0, 1.3);
     TH1D *jetPTresponse45 = new TH1D("jetPTreponse45", "jetPTresponse45", 50, 0, 1.3);
     TH1D *jetPTresponse46 = new TH1D("jetPTreponse46", "jetPTresponse46", 50, 0, 1.3);
     TH1D *jetPTresponse47 = new TH1D("jetPTreponse47", "jetPTresponse47", 50, 0, 1.3);
     TH1D *jetPTresponse48 = new TH1D("jetPTreponse48", "jetPTresponse48", 50, 0, 1.3);
     TH1D *jetPTresponse49 = new TH1D("jetPTreponse49", "jetPTresponse49", 50, 0, 1.3);


     TH1D *jetPTresponse50 = new TH1D("jetPTreponse50", "jetPTresponse50", 50, 0, 1.3);
     TH1D *jetPTresponse51 = new TH1D("jetPTreponse51", "jetPTresponse51", 50, 0, 1.3);
     TH1D *jetPTresponse52 = new TH1D("jetPTreponse52", "jetPTresponse52", 50, 0, 1.3);
     TH1D *jetPTresponse53 = new TH1D("jetPTreponse53", "jetPTresponse53", 50, 0, 1.3);
     TH1D *jetPTresponse54 = new TH1D("jetPTreponse54", "jetPTresponse54", 50, 0, 1.3);
     TH1D *jetPTresponse55 = new TH1D("jetPTreponse55", "jetPTresponse55", 50, 0, 1.3);
     TH1D *jetPTresponse56 = new TH1D("jetPTreponse56", "jetPTresponse56", 50, 0, 1.3);
     TH1D *jetPTresponse57 = new TH1D("jetPTreponse57", "jetPTresponse57", 50, 0, 1.3);
     TH1D *jetPTresponse58 = new TH1D("jetPTreponse58", "jetPTresponse58", 50, 0, 1.3);
     TH1D *jetPTresponse59 = new TH1D("jetPTreponse59", "jetPTresponse59", 50, 0, 1.3);


     TH1D *jetPTresponse60 = new TH1D("jetPTreponse60", "jetPTresponse60", 50, 0, 1.3);
     TH1D *jetPTresponse61 = new TH1D("jetPTreponse61", "jetPTresponse61", 50, 0, 1.3);
     TH1D *jetPTresponse62 = new TH1D("jetPTreponse62", "jetPTresponse62", 50, 0, 1.3);
     TH1D *jetPTresponse63 = new TH1D("jetPTreponse63", "jetPTresponse63", 50, 0, 1.3);
     TH1D *jetPTresponse64 = new TH1D("jetPTreponse64", "jetPTresponse64", 50, 0, 1.3);
     TH1D *jetPTresponse65 = new TH1D("jetPTreponse65", "jetPTresponse65", 50, 0, 1.3);
     TH1D *jetPTresponse66 = new TH1D("jetPTreponse66", "jetPTresponse66", 50, 0, 1.3);
     TH1D *jetPTresponse67 = new TH1D("jetPTreponse67", "jetPTresponse67", 50, 0, 1.3);
     TH1D *jetPTresponse68 = new TH1D("jetPTreponse68", "jetPTresponse68", 50, 0, 1.3);
     TH1D *jetPTresponse69 = new TH1D("jetPTreponse69", "jetPTresponse69", 50, 0, 1.3);


     TH1D *jetPTresponse70 = new TH1D("jetPTreponse70", "jetPTresponse70", 50, 0, 1.3);
     TH1D *jetPTresponse71 = new TH1D("jetPTreponse71", "jetPTresponse71", 50, 0, 1.3);
     TH1D *jetPTresponse72 = new TH1D("jetPTreponse72", "jetPTresponse72", 50, 0, 1.3);
     TH1D *jetPTresponse73 = new TH1D("jetPTreponse73", "jetPTresponse73", 50, 0, 1.3);
     TH1D *jetPTresponse74 = new TH1D("jetPTreponse74", "jetPTresponse74", 50, 0, 1.3);
     TH1D *jetPTresponse75 = new TH1D("jetPTreponse75", "jetPTresponse75", 50, 0, 1.3);
     TH1D *jetPTresponse76 = new TH1D("jetPTreponse76", "jetPTresponse76", 50, 0, 1.3);
     TH1D *jetPTresponse77 = new TH1D("jetPTreponse77", "jetPTresponse77", 50, 0, 1.3);
     TH1D *jetPTresponse78 = new TH1D("jetPTreponse78", "jetPTresponse78", 50, 0, 1.3);
     TH1D *jetPTresponse79 = new TH1D("jetPTreponse79", "jetPTresponse79", 50, 0, 1.3);


     TH1D *jetPTresponse80 = new TH1D("jetPTreponse80", "jetPTresponse80", 50, 0, 1.3);
     TH1D *jetPTresponse81 = new TH1D("jetPTreponse81", "jetPTresponse81", 50, 0, 1.3);
     TH1D *jetPTresponse82 = new TH1D("jetPTreponse82", "jetPTresponse82", 50, 0, 1.3);
     TH1D *jetPTresponse83 = new TH1D("jetPTreponse83", "jetPTresponse83", 50, 0, 1.3);
     TH1D *jetPTresponse84 = new TH1D("jetPTreponse84", "jetPTresponse84", 50, 0, 1.3);
     TH1D *jetPTresponse85 = new TH1D("jetPTreponse85", "jetPTresponse85", 50, 0, 1.3);
     TH1D *jetPTresponse86 = new TH1D("jetPTreponse86", "jetPTresponse86", 50, 0, 1.3);
     TH1D *jetPTresponse87 = new TH1D("jetPTreponse87", "jetPTresponse87", 50, 0, 1.3);
     TH1D *jetPTresponse88 = new TH1D("jetPTreponse88", "jetPTresponse88", 50, 0, 1.3);
     TH1D *jetPTresponse89 = new TH1D("jetPTreponse89", "jetPTresponse89", 50, 0, 1.3);


     TH1D *jetPTresponse90 = new TH1D("jetPTreponse90", "jetPTresponse90", 50, 0, 1.3);
     TH1D *jetPTresponse91 = new TH1D("jetPTreponse91", "jetPTresponse91", 50, 0, 1.3);
     TH1D *jetPTresponse92 = new TH1D("jetPTreponse92", "jetPTresponse92", 50, 0, 1.3);
     TH1D *jetPTresponse93 = new TH1D("jetPTreponse93", "jetPTresponse93", 50, 0, 1.3);
     TH1D *jetPTresponse94 = new TH1D("jetPTreponse94", "jetPTresponse94", 50, 0, 1.3);
     TH1D *jetPTresponse95 = new TH1D("jetPTreponse95", "jetPTresponse95", 50, 0, 1.3);
     TH1D *jetPTresponse96 = new TH1D("jetPTreponse96", "jetPTresponse96", 50, 0, 1.3);
     TH1D *jetPTresponse97 = new TH1D("jetPTreponse97", "jetPTresponse97", 50, 0, 1.3);
     TH1D *jetPTresponse98 = new TH1D("jetPTreponse98", "jetPTresponse98", 50, 0, 1.3);
     TH1D *jetPTresponse99 = new TH1D("jetPTreponse99", "jetPTresponse99", 50, 0, 1.3);

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
     Double_t AKTjetTheta;
     Double_t GenJetEta;
     Double_t GenJetPt;
     Double_t GenJetPhi;
     Double_t GenJetTheta;
     Double_t JER[10][10];
     Int_t    matchJetCnt[10][10];
     for (Int_t ThetaEntry=0; ThetaEntry < 10; ThetaEntry++) {
         for (Int_t PTentry=0; PTentry < 10; PTentry++) {
	     matchJetCnt[ThetaEntry][PTentry] = 0;
	     JER[ThetaEntry][PTentry]=0;
         } 
     }
     TH1D* jetPT2darray[10][10];

     jetPT2darray[0][0] = jetPTresponse00;
     jetPT2darray[0][1] = jetPTresponse01;
     jetPT2darray[0][2] = jetPTresponse02;
     jetPT2darray[0][3] = jetPTresponse03;
     jetPT2darray[0][4] = jetPTresponse04;
     jetPT2darray[0][5] = jetPTresponse05;
     jetPT2darray[0][6] = jetPTresponse06;
     jetPT2darray[0][7] = jetPTresponse07;
     jetPT2darray[0][8] = jetPTresponse08;
     jetPT2darray[0][9] = jetPTresponse09;

     jetPT2darray[1][0] = jetPTresponse10;
     jetPT2darray[1][1] = jetPTresponse11;
     jetPT2darray[1][2] = jetPTresponse12;
     jetPT2darray[1][3] = jetPTresponse13;
     jetPT2darray[1][4] = jetPTresponse14;
     jetPT2darray[1][5] = jetPTresponse15;
     jetPT2darray[1][6] = jetPTresponse16;
     jetPT2darray[1][7] = jetPTresponse17;
     jetPT2darray[1][8] = jetPTresponse18;
     jetPT2darray[1][9] = jetPTresponse19;


     jetPT2darray[2][0] = jetPTresponse20;
     jetPT2darray[2][1] = jetPTresponse21;
     jetPT2darray[2][2] = jetPTresponse22;
     jetPT2darray[2][3] = jetPTresponse23;
     jetPT2darray[2][4] = jetPTresponse24;
     jetPT2darray[2][5] = jetPTresponse25;
     jetPT2darray[2][6] = jetPTresponse26;
     jetPT2darray[2][7] = jetPTresponse27;
     jetPT2darray[2][8] = jetPTresponse28;
     jetPT2darray[2][9] = jetPTresponse29;


     jetPT2darray[3][0] = jetPTresponse30;
     jetPT2darray[3][1] = jetPTresponse31;
     jetPT2darray[3][2] = jetPTresponse32;
     jetPT2darray[3][3] = jetPTresponse33;
     jetPT2darray[3][4] = jetPTresponse34;
     jetPT2darray[3][5] = jetPTresponse35;
     jetPT2darray[3][6] = jetPTresponse36;
     jetPT2darray[3][7] = jetPTresponse37;
     jetPT2darray[3][8] = jetPTresponse38;
     jetPT2darray[3][9] = jetPTresponse39;


     jetPT2darray[4][0] = jetPTresponse40;
     jetPT2darray[4][1] = jetPTresponse41;
     jetPT2darray[4][2] = jetPTresponse42;
     jetPT2darray[4][3] = jetPTresponse43;
     jetPT2darray[4][4] = jetPTresponse44;
     jetPT2darray[4][5] = jetPTresponse45;
     jetPT2darray[4][6] = jetPTresponse46;
     jetPT2darray[4][7] = jetPTresponse47;
     jetPT2darray[4][8] = jetPTresponse48;
     jetPT2darray[4][9] = jetPTresponse49;


     jetPT2darray[5][0] = jetPTresponse50;
     jetPT2darray[5][1] = jetPTresponse51;
     jetPT2darray[5][2] = jetPTresponse52;
     jetPT2darray[5][3] = jetPTresponse53;
     jetPT2darray[5][4] = jetPTresponse54;
     jetPT2darray[5][5] = jetPTresponse55;
     jetPT2darray[5][6] = jetPTresponse56;
     jetPT2darray[5][7] = jetPTresponse57;
     jetPT2darray[5][8] = jetPTresponse58;
     jetPT2darray[5][9] = jetPTresponse59;


     jetPT2darray[6][0] = jetPTresponse60;
     jetPT2darray[6][1] = jetPTresponse61;
     jetPT2darray[6][2] = jetPTresponse62;
     jetPT2darray[6][3] = jetPTresponse63;
     jetPT2darray[6][4] = jetPTresponse64;
     jetPT2darray[6][5] = jetPTresponse65;
     jetPT2darray[6][6] = jetPTresponse66;
     jetPT2darray[6][7] = jetPTresponse67;
     jetPT2darray[6][8] = jetPTresponse68;
     jetPT2darray[6][9] = jetPTresponse69;


     jetPT2darray[7][0] = jetPTresponse70;
     jetPT2darray[7][1] = jetPTresponse71;
     jetPT2darray[7][2] = jetPTresponse72;
     jetPT2darray[7][3] = jetPTresponse73;
     jetPT2darray[7][4] = jetPTresponse74;
     jetPT2darray[7][5] = jetPTresponse75;
     jetPT2darray[7][6] = jetPTresponse76;
     jetPT2darray[7][7] = jetPTresponse77;
     jetPT2darray[7][8] = jetPTresponse78;
     jetPT2darray[7][9] = jetPTresponse79;


     jetPT2darray[8][0] = jetPTresponse80;
     jetPT2darray[8][1] = jetPTresponse81;
     jetPT2darray[8][2] = jetPTresponse82;
     jetPT2darray[8][3] = jetPTresponse83;
     jetPT2darray[8][4] = jetPTresponse84;
     jetPT2darray[8][5] = jetPTresponse85;
     jetPT2darray[8][6] = jetPTresponse86;
     jetPT2darray[8][7] = jetPTresponse87;
     jetPT2darray[8][8] = jetPTresponse88;
     jetPT2darray[8][9] = jetPTresponse89;


     jetPT2darray[9][0] = jetPTresponse90;
     jetPT2darray[9][1] = jetPTresponse91;
     jetPT2darray[9][2] = jetPTresponse92;
     jetPT2darray[9][3] = jetPTresponse93;
     jetPT2darray[9][4] = jetPTresponse94;
     jetPT2darray[9][5] = jetPTresponse95;
     jetPT2darray[9][6] = jetPTresponse96;
     jetPT2darray[9][7] = jetPTresponse97;
     jetPT2darray[9][8] = jetPTresponse98;
     jetPT2darray[9][9] = jetPTresponse99;


     jetPT2darray[0][0] = jetPTresponse00;
     jetPT2darray[0][1] = jetPTresponse01;
     jetPT2darray[0][2] = jetPTresponse02;
     jetPT2darray[0][3] = jetPTresponse03;
     jetPT2darray[0][4] = jetPTresponse04;
     jetPT2darray[0][5] = jetPTresponse05;
     jetPT2darray[0][6] = jetPTresponse06;
     jetPT2darray[0][7] = jetPTresponse07;
     jetPT2darray[0][8] = jetPTresponse08;
     jetPT2darray[0][9] = jetPTresponse09;


     Double_t JERtmp;
     Double_t DeltaPhi;
     Double_t DeltaR;
     Int_t PTGrid;
     Int_t ThetaGrid;
     cout << "Running Pairing up Algo..." << endl;
//Calculate JER
     for(Long64_t entry=0; entry < nEntries; entry++){
	 tree_sig->GetEntry(entry);
         AKTjet_size->GetBranch()->GetEntry(entry);
	 GenJet_size->GetBranch()->GetEntry(entry);

	 Int_t nAKTjet = AKTjet_size->GetValue();
	 Int_t nGenJet = GenJet_size->GetValue();

	 AKTjet_eta->GetBranch()->GetEntry(entry);
	 AKTjet_phi->GetBranch()->GetEntry(entry);
	 AKTjet_pt->GetBranch()->GetEntry(entry);
	 AKTjet_mass->GetBranch()->GetEntry(entry);

	 GenJet_eta->GetBranch()->GetEntry(entry);
	 GenJet_pt->GetBranch()->GetEntry(entry);
	 GenJet_mass->GetBranch()->GetEntry(entry);

	 TLorentzVector a;
	 TLorentzVector a1;
	 TLorentzVector a2;
	 a1.SetPtEtaPhiM(1,2,3,4);
	 a2.SetPtEtaPhiM(1,3,4,4);
         a = a1+a2;
         cout << a.Pt() << endl;
	 if (nAKTjet >= 4) {
             //Truth matching for all and calculate Jet energy response
	     for (Int_t aktentry=0; aktentry < nAKTjet; aktentry++) {
                 AKTjetEta = AKTjet_eta->GetValue(aktentry);
                 AKTjetPhi = AKTjet_phi->GetValue(aktentry);
	         AKTjetPt = AKTjet_pt->GetValue(aktentry);
	         AKTjetTheta = 2 * atan(exp(-AKTjetEta));
		 DeltaR = 100;
		 for (Int_t genentry=0; genentry < nGenJet; genentry++) {
		     GenJetEta = GenJet_eta->GetValue(genentry);
                     GenJetPhi = GenJet_phi->GetValue(genentry);
	             GenJetPt = GenJet_pt->GetValue(genentry);
	             GenJetTheta = 2 * atan(exp(-GenJetEta));
		     DeltaPhi = TMath::Abs(AKTjetPhi-GenJetPhi);
		     if (DeltaPhi > TMath::Pi()) {
		         DeltaPhi -= TMath::TwoPi();
		     }
		     Float_t DeltaRtmp = TMath::Sqrt(pow((AKTjetEta-GenJetEta),2)+pow(DeltaPhi,2));
		     if (DeltaRtmp < DeltaR) {
		        DeltaR = DeltaRtmp;
			JERtmp = AKTjetPt/GenJetPt;
		     }
		 }

		 if (DeltaR < 0.5) {
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
		     if (1.8< AKTjetTheta and AKTjetTheta <= 2.1) {
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
		     matchJetCnt[ThetaGrid][PTGrid]+=1;
                     JER[ThetaGrid][PTGrid] = (JER[ThetaGrid][PTGrid] * (matchJetCnt[ThetaGrid][PTGrid]-1)+JERtmp)/matchJetCnt[ThetaGrid][PTGrid];
		     jetPT2darray[ThetaGrid][PTGrid]->Fill(JERtmp);
		 } else {
		     //cout << DeltaR << endl;
		 }
	     }
	 }
     }
     TCanvas *mycanvas = new TCanvas("mycanvas","My Canvas",800,600);
//Calibration result
     cout << "Jet Energy Calibration gives following Jet Energy Response:" << endl;
     for (Int_t PTentry=9; PTentry >=0; PTentry-=1) {
         for (Int_t ThetaEntry=0; ThetaEntry < 10; ThetaEntry++) {
	     if (matchJetCnt[ThetaEntry][PTentry] == 0){
	         JER[ThetaEntry][PTentry] = 1.0000;
	     }
             jetPT2darray[ThetaEntry][PTentry]->GetXaxis()->SetTitle("Jet P_{T} Response");
             jetPT2darray[ThetaEntry][PTentry]->Draw();
	     mycanvas->SaveAs(TString::Format("jetPTresponse_%d%d.png",ThetaEntry,PTentry));
             jetPT2darray[ThetaEntry][PTentry]->Write();
         } 
     }
     for (Int_t PTentry=9; PTentry >=0; PTentry-=1) {
         for (Int_t ThetaEntry=0; ThetaEntry < 10; ThetaEntry++) {
	     cout << floorf(JER[ThetaEntry][PTentry]*10000)/10000<< ", ";
         } 
	 cout << endl;
     }

//run over event entries
     for(Long64_t entry=0; entry < nEntries; entry++){
	 tree_sig->GetEntry(entry);
	 tree_output->GetEntry(entry);
         AKTjet_size->GetBranch()->GetEntry(entry);
	 GenJet_size->GetBranch()->GetEntry(entry);

	 Int_t nAKTjet = AKTjet_size->GetValue();
	 Int_t nGenJet = GenJet_size->GetValue();

	 AKTjet_eta->GetBranch()->GetEntry(entry);
	 AKTjet_phi->GetBranch()->GetEntry(entry);
	 AKTjet_pt->GetBranch()->GetEntry(entry);
	 AKTjet_mass->GetBranch()->GetEntry(entry);

	 GenJet_eta->GetBranch()->GetEntry(entry);
	 GenJet_phi->GetBranch()->GetEntry(entry);
	 GenJet_pt->GetBranch()->GetEntry(entry);
	 GenJet_mass->GetBranch()->GetEntry(entry);
         
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

         Double_t AKTjetpair[nAKTjet*(nAKTjet-1)/2][6];

         AKTjetpaircnt = 0;
	 
	 if (nAKTjet >= 4) {
             for (Int_t aktentry=0; aktentry < nAKTjet; aktentry++) {
                 AKTjetEta = AKTjet_eta->GetValue(aktentry);
	         AKTjetPt = AKTjet_pt->GetValue(aktentry);
                 AKTjetTheta = 2 * atan(exp(-AKTjetEta));
	         AKTjetPT_Theta->Fill(AKTjetTheta, AKTjetPt);
		 JetEnergyFix(AKTjetEta, AKTjetPt, JER);
	     }
             for (Int_t genentry=0; genentry < nGenJet; genentry++) {
                 GenJetEta = GenJet_eta->GetValue(genentry);
	         GenJetPt = GenJet_pt->GetValue(genentry);
	         GenJetTheta = 2 * atan(exp(-GenJetEta));
	         GenJetPT_Theta->Fill(GenJetTheta, GenJetPt);
	     }
             //Pairing up leading anti-KT jet pair
	     for (Int_t akt1entry=0; akt1entry < nAKTjet; akt1entry++){
	         for (Int_t akt2entry=akt1entry+1; akt2entry < nAKTjet; akt2entry++){

		     AKTjet1eta1 = AKTjet_eta->GetValue(akt1entry);
         	     AKTjet1phi1 = AKTjet_phi->GetValue(akt1entry);
         	     AKTjet1pt1 = AKTjet_pt->GetValue(akt1entry);
		     JetEnergyFix(AKTjet1eta1, AKTjet1pt1, JER);
         	     AKTjet1mass1 = AKTjet_mass->GetValue(akt1entry);
                     AKTjet1eta2 = AKTjet_eta->GetValue(akt2entry);
         	     AKTjet1phi2 = AKTjet_phi->GetValue(akt2entry);
         	     AKTjet1pt2 = AKTjet_pt->GetValue(akt2entry);
		     JetEnergyFix(AKTjet1eta2, AKTjet1pt2, JER);
         	     AKTjet1mass2 = AKTjet_mass->GetValue(akt2entry);
		     AKTjetpairmass = 0;

		     AKTjet1.SetPtEtaPhiM(AKTjet1pt1, AKTjet1eta1, AKTjet1phi1,AKTjet1mass1);
		     AKTjet2.SetPtEtaPhiM(AKTjet1pt2, AKTjet1eta2, AKTjet1phi2,AKTjet1mass2);
		     AKTh1=AKTjet1+AKTjet2;
		     AKTjetpairmass = AKTh1.Mag();

		     AKTjetpair[AKTjetpaircnt][0]=AKTjetpairmass;
		     AKTjetpair[AKTjetpaircnt][1]=akt1entry;
		     AKTjetpair[AKTjetpaircnt][2]=akt2entry;
                     
		     Int_t akt1entry2;
		     Int_t akt2entry2;
		     for (akt1entry2 = 0; akt1entry2 < nAKTjet; akt1entry2++){
	                 if ((akt1entry2 != akt1entry) and (akt1entry2 != akt2entry)){
		             for (akt2entry2 = 0; akt2entry2 < nAKTjet; akt2entry2++){
		                 if (((akt2entry2 != akt1entry) and (akt2entry2 != akt2entry)) and (akt2entry2 != akt1entry)){
			      
			             AKTjet2eta1 = AKTjet_eta->GetValue(akt1entry2);
         	                     AKTjet2phi1 = AKTjet_phi->GetValue(akt1entry2);
         	                     AKTjet2pt1 = AKTjet_pt->GetValue(akt1entry2);
		                     JetEnergyFix(AKTjet2eta1, AKTjet2pt1, JER);
         	                     AKTjet2mass1 = AKTjet_mass->GetValue(akt1entry2);
                                     AKTjet2eta2 = AKTjet_eta->GetValue(akt2entry2);
         	                     AKTjet2phi2 = AKTjet_phi->GetValue(akt2entry2);
         	                     AKTjet2pt2 = AKTjet_pt->GetValue(akt2entry2);
		                     JetEnergyFix(AKTjet2eta2, AKTjet2pt2, JER);
         	                     AKTjet2mass2 = AKTjet_mass->GetValue(akt2entry2);
		                     AKTjetpairmass = 0;

		                     AKTjet1.SetPtEtaPhiM(AKTjet2pt1, AKTjet2eta1, AKTjet2phi1, AKTjet2mass1);
		                     AKTjet2.SetPtEtaPhiM(AKTjet2pt2, AKTjet2eta2, AKTjet2phi2, AKTjet2mass2);
		                     AKTh2=AKTjet1+AKTjet2;
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
	             
                     AKTjetpair[AKTjetpaircnt][3]=AKTjetpair2Mass;
		     AKTjetpair[AKTjetpaircnt][4]=AKTjet2entry1;
		     AKTjetpair[AKTjetpaircnt][5]=AKTjet2entry2;
	         }  

		 AKTjetpaircnt += 1;
	     }
             //cout << AKTjetpaircnt;
	     for (Int_t entry = 0; entry < (nAKTjet-1)*nAKTjet/2; entry++) {
	         if (abs(125 - AKTjetpair[entry][0])*abs(125 - AKTjetpair[entry][3]) < diHiggsdis) {
                     
		     pair1jet1entry =  AKTjetpair[entry][1];
                     pair1jet2entry =  AKTjetpair[entry][2];
                     pair2jet1entry =  AKTjetpair[entry][4];
                     pair2jet2entry =  AKTjetpair[entry][5];

		     AKTjetpair1Mass = AKTjetpair[entry][0];
		     AKTjetpair2Mass = AKTjetpair[entry][3];
		     
		     diHiggsdis = abs(125 - AKTjetpair[entry][0])*abs(125 - AKTjetpair[entry][3]);
		 }		 
	     }
	      
	     AKTjet1eta1 = AKTjet_eta->GetValue(pair1jet1entry);
             AKTjet1phi1 = AKTjet_phi->GetValue(pair1jet1entry);
             AKTjet1pt1 = AKTjet_pt->GetValue(pair1jet1entry);
             JetEnergyFix(AKTjet1eta1, AKTjet1pt1, JER);
             AKTjet1mass1 = AKTjet_mass->GetValue(pair1jet1entry);
             AKTjet1eta2 = AKTjet_eta->GetValue(pair1jet2entry);
             AKTjet1phi2 = AKTjet_phi->GetValue(pair1jet2entry);
             AKTjet1pt2 = AKTjet_pt->GetValue(pair1jet2entry);
             JetEnergyFix(AKTjet1eta2, AKTjet1pt2, JER);
             AKTjet1mass2 = AKTjet_mass->GetValue(pair1jet2entry);

             AKTjet2eta1 = AKTjet_eta->GetValue(pair2jet1entry);
             AKTjet2phi1 = AKTjet_phi->GetValue(pair2jet1entry);
             AKTjet2pt1 = AKTjet_pt->GetValue(pair2jet1entry);
             JetEnergyFix(AKTjet2eta1, AKTjet2pt1, JER);
             AKTjet2mass1 = AKTjet_mass->GetValue(pair2jet1entry);
             AKTjet2eta2 = AKTjet_eta->GetValue(pair2jet2entry);
             AKTjet2phi2 = AKTjet_phi->GetValue(pair2jet2entry);
             AKTjet2pt2 = AKTjet_pt->GetValue(pair2jet2entry);
             JetEnergyFix(AKTjet2eta2, AKTjet2pt2, JER);
             AKTjet2mass2 = AKTjet_mass->GetValue(pair2jet2entry);

	     bool AKT1jet1flag = false;
	     bool AKT1jet2flag = false;
             bool AKT2jet1flag = false;
	     bool AKT2jet2flag = false;

	     bool AKT1jet1edgeflag = false;
	     bool AKT1jet2edgeflag = false;
	     bool AKT2jet1edgeflag = false;
	     bool AKT2jet2edgeflag = false;

             jet1DeltaR = 100;
             jet2DeltaR = 100;
             jet3DeltaR = 100;
             jet4DeltaR = 100;

             Float_t jet1DeltaPhi;
             Float_t jet2DeltaPhi;
             Float_t jet3DeltaPhi;
             Float_t jet4DeltaPhi;

	     Int_t jet1entry;
	     Int_t jet2entry;
	     Int_t jet3entry;
	     Int_t jet4entry;
//Truth matching for anti-KT jets
             for (Int_t gen1entry=0; gen1entry < nGenJet; gen1entry++){
  	         Gen1eta = GenJet_eta->GetValue(gen1entry);
	         Gen1phi = GenJet_phi->GetValue(gen1entry);
		 jet1DeltaPhi = std::abs(AKTjet1phi1-Gen1phi);                
		 if (jet1DeltaPhi > TMath::Pi()) {
		     jet1DeltaPhi -= 2 * TMath::Pi();
		 } 
	         Float_t jet1DeltaRtmp = TMath::Sqrt(pow((AKTjet1eta1-Gen1eta),2)+pow(jet1DeltaPhi,2));
	         if (jet1DeltaRtmp < jet1DeltaR){
	             jet1DeltaR = jet1DeltaRtmp;
		     if (jet1DeltaR < 0.5){
		         //if (abs(Gen1eta) < 2.25){
		             AKT1jet1flag = true;
		             jet1entry = gen1entry;
		         //}
		     }
	         }  
             }
             for (Int_t gen2entry=0; gen2entry < nGenJet; gen2entry++){
  	         Gen2eta = GenJet_eta->GetValue(gen2entry);
	         Gen2phi = GenJet_phi->GetValue(gen2entry);
                 jet2DeltaPhi = std::abs(AKTjet1phi2-Gen2phi);                
		 if (jet2DeltaPhi > TMath::Pi()) {
		     jet2DeltaPhi -= 2 * TMath::Pi();
		 } 
	         Float_t jet2DeltaRtmp = TMath::Sqrt(pow((AKTjet1eta2-Gen2eta),2)+pow(jet2DeltaPhi,2));
	         if (gen2entry != jet1entry){
                     if (jet2DeltaRtmp < jet2DeltaR){
	                 jet2DeltaR = jet2DeltaRtmp;
                         if (jet2DeltaR < 0.5){
		             //if (abs(Gen2eta) < 2.25){
		                 AKT1jet2flag = true;
		                 jet2entry = gen2entry;
			     //}
		         }
	             }
	         }
             }

	     for (Int_t gen3entry=0; gen3entry < nGenJet; gen3entry++){
  	         Gen3eta = GenJet_eta->GetValue(gen3entry);
	         Gen3phi = GenJet_phi->GetValue(gen3entry);
                 jet3DeltaPhi = std::abs(AKTjet2phi1-Gen3phi);                
		 if (jet3DeltaPhi > TMath::Pi()) {
		     jet3DeltaPhi -= 2 * TMath::Pi();
		 } 
	         Float_t jet1DeltaRtmp = TMath::Sqrt(pow((AKTjet2eta1-Gen3eta),2)+pow(jet3DeltaPhi,2));
	         if (jet1DeltaRtmp < jet3DeltaR){
	             jet3DeltaR = jet1DeltaRtmp;
		     if (jet3DeltaR < 0.5){
		         //if (abs(Gen1eta) < 2.25){
		             AKT2jet1flag = true;
		             jet3entry = gen3entry;
		         //}
		     }
	         }  
             }
             for (Int_t gen4entry=0; gen4entry < nGenJet; gen4entry++){
  	         Gen4eta = GenJet_eta->GetValue(gen4entry);
	         Gen4phi = GenJet_phi->GetValue(gen4entry);
                 jet4DeltaPhi = std::abs(AKTjet2phi2-Gen4phi);                
		 if (jet4DeltaPhi > TMath::Pi()) {
		     jet4DeltaPhi -= 2 * TMath::Pi();
		 } 
	         Float_t jet2DeltaRtmp = TMath::Sqrt(pow((AKTjet2eta2-Gen4eta),2)+pow(jet4DeltaPhi,2));
	         if (gen4entry != jet1entry){
                     if (jet2DeltaRtmp < jet4DeltaR){
	                 jet4DeltaR = jet2DeltaRtmp;
                         if (jet4DeltaR < 0.5){
		             //if (abs(Gen2eta) < 2.25){
		                 AKT2jet2flag = true;
		                 jet4entry = gen4entry;
			     //}
		         }
	             }
	         }
             }
	     
	     if (AKTjetpair1Mass < AKTjetpair2Mass) {
	         swap(AKTjetpair1Mass, AKTjetpair2Mass);
		 swap(jet1entry,jet3entry);
		 swap(jet2entry,jet4entry);
		 swap(AKTjet1pt1,AKTjet2pt1);
		 swap(AKTjet2pt1,AKTjet2pt2);
		 swap(AKTjet1eta1,AKTjet2eta1);
		 swap(AKTjet1eta2,AKTjet2eta2);
	     }
	      
             Gen1eta = GenJet_eta->GetValue(jet1entry);
             Gen1phi = GenJet_phi->GetValue(jet1entry);
             Gen1pt = GenJet_pt->GetValue(jet1entry);
             Gen1mass = GenJet_mass->GetValue(jet1entry);
             Gen2eta = GenJet_eta->GetValue(jet2entry);
             Gen2phi = GenJet_phi->GetValue(jet2entry);
             Gen2pt = GenJet_pt->GetValue(jet2entry);
             Gen2mass = GenJet_mass->GetValue(jet2entry);
         
	     GenJet1.SetPtEtaPhiM(Gen1pt, Gen1eta, Gen1phi,Gen1mass);
	     GenJet2.SetPtEtaPhiM(Gen2pt, Gen2eta, Gen2phi,Gen2mass);
	     Genh2=GenJet1+GenJet2;
	     GenJetMass = Genh2.Mag();
             Double_t AKTGenMass1diff = (AKTjetpair1Mass - GenJetMass)/GenJetMass;

	     Gen3eta = GenJet_eta->GetValue(jet3entry);
             Gen3phi = GenJet_phi->GetValue(jet3entry);
             Gen3pt = GenJet_pt->GetValue(jet3entry);
             Gen3mass = GenJet_mass->GetValue(jet3entry);
             Gen4eta = GenJet_eta->GetValue(jet4entry);
             Gen4phi = GenJet_phi->GetValue(jet4entry);
             Gen4pt = GenJet_pt->GetValue(jet4entry);
             Gen4mass = GenJet_mass->GetValue(jet4entry);
         
	     GenJet3.SetPtEtaPhiM(Gen3pt, Gen3eta, Gen3phi,Gen3mass);
	     GenJet4.SetPtEtaPhiM(Gen4pt, Gen4eta, Gen4phi,Gen4mass);
	     Genh2=GenJet3+GenJet4;
	     GenJetMass = Genh2.Mag();
             Double_t AKTGenMass2diff = (AKTjetpair2Mass - GenJetMass)/GenJetMass;
             
	     Double_t AKTpt[] = { AKTjet1pt1, AKTjet1pt2, AKTjet2pt1, AKTjet2pt2 };
	     //int AKTsize = sizeof(AKTpt) / sizeof(AKTpt[0]);
             
	     Double_t Genpt[] = { Gen1pt, Gen2pt, Gen3pt, Gen4pt };
	     //int Gensize = sizeof(Genpt) / sizeof(Genpt[0]);

	     Double_t deltaR[] = { jet1DeltaR, jet2DeltaR, jet3DeltaR, jet4DeltaR };

	     Double_t AKTeta[] = { AKTjet1eta1, AKTjet1eta2, AKTjet2eta1, AKTjet2eta2 };
             
	     //sort pt 
	     for (int i=0; i<3; i++) {
	         for (int j=i+1; j<4; j++) {
		     if (AKTpt[j] > AKTpt[i]) {
                         swap(AKTpt[j], AKTpt[i]);
			 swap(Genpt[j], Genpt[i]);
			 swap(deltaR[j], deltaR[i]);
			 swap(AKTeta[j], AKTeta[i]);
		     }
		 }
	     }
	     //calculate resolution
	     Double_t AKTGenPt1diff = (AKTpt[0] - Genpt[0])/Genpt[0];
	     Double_t AKTGenPt2diff = (AKTpt[1] - Genpt[1])/Genpt[1];
	     Double_t AKTGenPt3diff = (AKTpt[2] - Genpt[2])/Genpt[2];
	     Double_t AKTGenPt4diff = (AKTpt[3] - Genpt[3])/Genpt[3];
	     
	     //calculate angle
             AKTjet1theta1 = 2 * atan(exp(-AKTeta[0]));
             AKTjet1theta2 = 2 * atan(exp(-AKTeta[1]));
             AKTjet2theta1 = 2 * atan(exp(-AKTeta[2]));
             AKTjet2theta2 = 2 * atan(exp(-AKTeta[3]));
             
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
	     
             //if (AKT1jet1flag==true and AKT1jet2flag==true and AKT2jet1flag==true and AKT2jet2flag==true){
	         AKTGenPt1Comp->Fill(AKTGenPt1diff);
		 AKTGenPt2Comp->Fill(AKTGenPt2diff);
		 AKTGenPt3Comp->Fill(AKTGenPt3diff);
		 AKTGenPt4Comp->Fill(AKTGenPt4diff);
		 alljet1DeltaR->Fill(deltaR[0]);
		 alljet2DeltaR->Fill(deltaR[1]);
		 alljet3DeltaR->Fill(deltaR[2]);
		 alljet4DeltaR->Fill(deltaR[3]);
		 alljet1theta->Fill(AKTjet1theta1);
		 alljet2theta->Fill(AKTjet1theta2);
		 alljet3theta->Fill(AKTjet2theta1);
		 alljet4theta->Fill(AKTjet2theta2);
		 /*
                 alljet1eta->Fill(AKTeta[0]);
		 alljet2eta->Fill(AKTeta[1]);
		 alljet3eta->Fill(AKTeta[2]);
		 alljet4eta->Fill(AKTeta[3]);
		 */
		 if (abs(AKTGenPt1diff)>0.2) {
		     jet1Reso_Pt->Fill(Genpt[0], AKTGenPt1diff);
                     jet1Reso_DeltaR->Fill(deltaR[0], AKTGenPt1diff);
		     badjet1DeltaR->Fill(deltaR[0]);
		     //badjet1eta->Fill(AKTeta[0]);
		     badjet1theta->Fill(AKTjet1theta1);
		 } else {
		     goodjet1theta->Fill(AKTjet1theta1);
		     //goodjet1eta->Fill(AKTeta[0]);
		 }
                 if (abs(AKTGenPt2diff)>0.2) {
		     jet2Reso_Pt->Fill(Genpt[1], AKTGenPt2diff);
		     jet2Reso_DeltaR->Fill(deltaR[1], AKTGenPt2diff);
		     badjet2DeltaR->Fill(deltaR[1]);
		     //badjet2eta->Fill(AKTeta[1]);
		     badjet2theta->Fill(AKTjet1theta2);
		 } else {
		     goodjet2theta->Fill(AKTjet1theta2);
		     //goodjet2eta->Fill(AKTeta[1]);
		 }
                 if (abs(AKTGenPt3diff)>0.2) {
		     jet3Reso_Pt->Fill(Genpt[2], AKTGenPt3diff);
		     jet3Reso_DeltaR->Fill(deltaR[2], AKTGenPt3diff);
		     badjet3DeltaR->Fill(deltaR[2]);
		     //badjet3eta->Fill(AKTeta[2]);
		     badjet3theta->Fill(AKTjet2theta1);
		 } else {
		     goodjet3theta->Fill(AKTjet2theta1);
		     //goodjet3eta->Fill(AKTeta[2]);
		 }
                 if (abs(AKTGenPt4diff)>0.2) {
		     jet4Reso_Pt->Fill(Genpt[3], AKTGenPt4diff);
		     jet4Reso_DeltaR->Fill(deltaR[3], AKTGenPt4diff);
		     badjet4DeltaR->Fill(deltaR[3]);
		     //badjet4eta->Fill(AKTeta[3]);
		     badjet4theta->Fill(AKTjet2theta2);
		 } else {
		     goodjet4theta->Fill(AKTjet2theta2);
		     //goodjet4eta->Fill(AKTeta[3]);
		 }


	         if (abs(AKTGenPt1diff)<0.15 and abs(AKTGenPt2diff)<0.15 and abs(AKTGenPt3diff)<0.15 and abs(AKTGenPt4diff)<0.15) { 
		     //if (AKT1jet1edgeflag==true and AKT1jet2edgeflag==true and AKT2jet1edgeflag==true and AKT2jet2edgeflag==true) {
		         AKTjetMass1->Fill(AKTjetpair1Mass);
                         AKTjetMass2->Fill(AKTjetpair2Mass);
		         //GenAKTMass2->Fill(GenJetMass);
	                 AKTGenMass1Comp->Fill(AKTGenMass1diff);
	                 AKTGenMass2Comp->Fill(AKTGenMass2diff);
		     //}
	         }
	     //} 
             //AKTjetMass1->Fill(AKTjetpair1Mass);
             //AKTjetMass2->Fill(AKTjetpair2Mass);

         }
	 //GenUncutMass2->Fill(GenJetMass);
     }

 //Fitting and plotting
     TF1 *jetpair1fit = new TF1("jetpair1fit", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)",25,600);
     TF1 *jetpair2fit = new TF1("jetpair2fit", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)+expo(5)",25,600);
     TF1 *fSignal = new TF1("fSignal","gaus+gaus(3)",20,600);
     TF1 *fBackground = new TF1("fBackground","expo",20,600);
     TF1 *f1peak = new TF1("f1peak","gaus",20,600);
     TF1 *f2peak = new TF1("f2peak","gaus",20,600);
     Double_t param[8];

     jetpair2fit->SetParameters(200,120,10,20,10,2,-0.0001);
     jetpair2fit->SetParLimits(0,30,100);
     jetpair2fit->SetParLimits(1,100,140);
     jetpair2fit->SetParLimits(2,10,30);
     jetpair2fit->SetParLimits(5,0,20);
     jetpair2fit->SetParLimits(6,-1,-0.0001);
     //jetpair2fit->SetParLimits(4,50,109);
     jetpair2fit->SetParLimits(3,10,50);
     jetpair2fit->SetParLimits(4,0,10);
     
     jetpair1fit->SetParameters(300,120,10,40,20);
     jetpair1fit->SetParLimits(0,50,350);
     jetpair1fit->SetParLimits(1,100,140);
     jetpair1fit->SetParLimits(2,5,30);
     /*jetpair1fit->SetParLimits(6,0,8);
     jetpair1fit->SetParLimits(7,-1.5,-0.0001);*/
     jetpair1fit->SetParLimits(3,0,60);
     jetpair1fit->SetParLimits(4,5,40);
	     
     cout <<endl<< "Run gaussian fit for the leading anti-KT jet pair..."<<endl;
     AKTjetMass1->Fit("jetpair1fit","R");
     jetpair1fit->GetParameters(param);
     fSignal->SetParameters(&param[0]);
     fBackground->SetParameters(&param[6]);
     f1peak->SetParameters(&param[0]);
     f2peak->SetParameters(param[3],param[1],param[4]);
     TH1D *AKTjetMass1Signal = new TH1D(*AKTjetMass1);
     AKTjetMass1Signal->Sumw2();
     AKTjetMass1Signal->Add(fBackground,-1);
     AKTjetMass1->GetXaxis()->SetTitle("M [GeV]");
     AKTjetMass1->GetYaxis()->SetTitle("Events");
     AKTjetMass1->Draw();fBackground->SetLineColor(4);fBackground->Draw("SAME");
     f1peak->SetLineColor(3);f1peak->Draw("SAME");
     f2peak->SetLineColor(6);f2peak->Draw("SAME");
     mycanvas->SaveAs("AKTjetpair1Mass.png");

     AKTjetMass2->Fit("jetpair2fit","R");
     jetpair2fit->GetParameters(param);
     fSignal->SetParameters(&param[0]);
     fBackground->SetParameters(&param[5]);
     f1peak->SetParameters(&param[0]);
     f2peak->SetParameters(param[3],param[1],param[4]);
     TH1D *AKTjetMass2Signal = new TH1D(*AKTjetMass2);
     AKTjetMass2Signal->Sumw2();
     AKTjetMass2Signal->Add(fBackground,-1);
     AKTjetMass2->GetXaxis()->SetTitle("M [GeV]");
     AKTjetMass2->GetYaxis()->SetTitle("Events");
     AKTjetMass2->Draw();fBackground->SetLineColor(4);fBackground->Draw("SAME"); 
     f1peak->SetLineColor(3);f1peak->Draw("SAME");
     f2peak->SetLineColor(6);f2peak->Draw("SAME");
     //AKTjetMass2Signal->Draw("SAME"); fSignal->Draw("SAME");
     mycanvas->SaveAs("AKTjetpair2Mass.png");
     cout <<endl<< "Run gaussian fit for the sub-leading GenJet pair under selection of AKT..."<<endl;
     //GenAKTMass2->Fit("jetpair2fit","R");
     //mycanvas->SaveAs("GenAKTMass2.png");
     cout <<endl<< "Run gaussian fit for the sub-leading GenJet pair (Uncut)..."<<endl;
     //GenUncutMass2->Fit("jetpair2fit","R");
     //mycanvas->SaveAs("GenUncutMass2.png");
     
     AKTGenMass1Comp->GetXaxis()->SetTitle("(M_{H1}-M_{Gen})/M_{Gen}");
     AKTGenMass1Comp->GetYaxis()->SetTitle("Events");
     AKTGenMass1Comp->Draw();
     mycanvas->SaveAs("AKTGenMass1Comp.png");
     
     AKTGenMass2Comp->GetXaxis()->SetTitle("(M_{H2}-M_{Gen})/M_{Gen}");
     AKTGenMass2Comp->GetYaxis()->SetTitle("Events");
     AKTGenMass2Comp->Draw();
     mycanvas->SaveAs("AKTGenMass2Comp.png");

     AKTGenPt1Comp->GetXaxis()->SetTitle("(P_{T_1}-P_{T_{Gen1}})/P_{T_{Gen1}}");
     AKTGenPt1Comp->GetYaxis()->SetTitle("Events");
     AKTGenPt1Comp->Draw();
     mycanvas->SaveAs("AKTGenPt1Comp.png");
 
     AKTGenPt2Comp->GetXaxis()->SetTitle("(P_{T_2}-P_{T_{Gen2}})/P_{T_{Gen2}}");
     AKTGenPt2Comp->GetYaxis()->SetTitle("Events");
     AKTGenPt2Comp->Draw();
     mycanvas->SaveAs("AKTGenPt2Comp.png");
 
     AKTGenPt3Comp->GetXaxis()->SetTitle("(P_{T_3}-P_{T_{Gen3}})/P_{T_{Gen3}}");
     AKTGenPt3Comp->GetYaxis()->SetTitle("Events");
     AKTGenPt3Comp->Draw();
     mycanvas->SaveAs("AKTGenPt3Comp.png");
 
     AKTGenPt4Comp->GetXaxis()->SetTitle("(P_{T_4}-P_{T_{Gen4}})/P_{T_{Gen4}}");
     AKTGenPt4Comp->GetYaxis()->SetTitle("Events");
     AKTGenPt4Comp->Draw();
     mycanvas->SaveAs("AKTGenPt4Comp.png");
      
     //Calibration plot
     GenJetPT_Theta->GetXaxis()->SetTitle("#theta");
     GenJetPT_Theta->GetYaxis()->SetTitle("P_{T} [GeV]");
     GenJetPT_Theta->SetMarkerColor(kBlue);
     //GenJetPT_Theta->SetMarkerStyle(kCircle);
     GenJetPT_Theta->SetLineColor(kBlue);
     GenJetPT_Theta->Draw();
     AKTjetPT_Theta->SetMarkerColor(kRed);
     //AKTjetPT_Theta->SetMarkerStyle(kPlus);
     AKTjetPT_Theta->SetLineColor(kRed);
     AKTjetPT_Theta->Draw("same");
     mycanvas->SaveAs("Calibration.png");

     jet1Reso_Pt->GetXaxis()->SetTitle("P_{T} [GeV]");
     jet1Reso_Pt->GetYaxis()->SetTitle("Resolution (P_{T})");
     jet1Reso_Pt->Draw();//"COLZ");
     mycanvas->SaveAs("jet1Reso_Pt.png");
     
     jet2Reso_Pt->GetXaxis()->SetTitle("P_{T} [GeV]");
     jet2Reso_Pt->GetYaxis()->SetTitle("Resolution (P_{T})");
     jet2Reso_Pt->Draw();//"COLZ");
     mycanvas->SaveAs("jet2Reso_Pt.png");

     jet3Reso_Pt->GetXaxis()->SetTitle("P_{T} [GeV]");
     jet3Reso_Pt->GetYaxis()->SetTitle("Resolution (P_{T})");
     jet3Reso_Pt->Draw();//"COLZ");
     mycanvas->SaveAs("jet3Reso_Pt.png");
 
     jet4Reso_Pt->GetXaxis()->SetTitle("P_{T} [GeV]");
     jet4Reso_Pt->GetYaxis()->SetTitle("Resolution (P_{T})");
     jet4Reso_Pt->Draw();//"COLZ");
     mycanvas->SaveAs("jet4Reso_Pt.png");
     
     jet1Reso_DeltaR->GetXaxis()->SetTitle("#Delta_{R}");
     jet1Reso_DeltaR->GetYaxis()->SetTitle("Resolution (P_{T})");   
     jet1Reso_DeltaR->Draw();//"COLZ");
     mycanvas->SaveAs("jet1Reso_DeltaR.png");
 
     jet2Reso_DeltaR->GetXaxis()->SetTitle("#Delta_{R}");
     jet2Reso_DeltaR->GetYaxis()->SetTitle("Resolution (P_{T})"); 
     jet2Reso_DeltaR->Draw();//"COLZ");
     mycanvas->SaveAs("jet2Reso_DeltaR.png");

     jet3Reso_DeltaR->GetXaxis()->SetTitle("#Delta_{R}");
     jet3Reso_DeltaR->GetYaxis()->SetTitle("Resolution (P_{T})"); 
     jet3Reso_DeltaR->Draw();//"COLZ");
     mycanvas->SaveAs("jet3Reso_DeltaR.png");

     jet4Reso_DeltaR->GetXaxis()->SetTitle("#Delta_{R}");
     jet4Reso_DeltaR->GetYaxis()->SetTitle("Resolution (P_{T})");     
     jet4Reso_DeltaR->Draw();//"COLZ");
     mycanvas->SaveAs("jet4Reso_DeltaR.png");
/*
     alljet1eta->Draw();
     alljet1eta->SetLineColor(kRed);
     alljet1eta->GetXaxis()->SetTitle("#eta");
     alljet1eta->GetYaxis()->SetTitle("Events");   
     badjet1eta->SetStats(0);
     badjet1eta->Draw("same");
     badjet1eta->SetLineColor(kBlue+1);
     TH1D *jet1etaRatio = (TH1D*)badjet1eta->Clone("jet1etaRatio");
     jet1etaRatio->SetLineColor(kBlack);
     jet1etaRatio->Sumw2();
     jet1etaRatio->SetStats(0);
     jet1etaRatio->Divide(alljet1eta);
     jet1etaRatio->Scale(100);
     goodjet1eta->Draw("same");
     goodjet1eta->SetLineColor(kGreen);
     mycanvas->SaveAs("AKTjet1eta.png");
     jet1etaRatio->Draw("HIST");
     mycanvas->SaveAs("jet1etaRatio.png");
 
     alljet2eta->Draw();
     alljet2eta->SetLineColor(kRed);
     alljet2eta->GetXaxis()->SetTitle("#eta");
     alljet2eta->GetYaxis()->SetTitle("Events");   
     badjet2eta->SetStats(0);
     badjet2eta->Draw("same");
     badjet2eta->SetLineColor(kBlue+1);
     TH1D *jet2etaRatio = (TH1D*)badjet2eta->Clone("jet2etaRatio");
     jet2etaRatio->SetLineColor(kBlack);
     jet2etaRatio->Sumw2();
     jet2etaRatio->SetStats(0);
     jet2etaRatio->Divide(alljet2eta);
     jet2etaRatio->Scale(100);
     goodjet2eta->Draw("same");
     goodjet2eta->SetLineColor(kGreen);
     mycanvas->SaveAs("AKTjet2eta.png");
     jet2etaRatio->Draw("HIST");
     mycanvas->SaveAs("jet2etaRatio.png");
 
     alljet3eta->Draw();
     alljet3eta->SetLineColor(kRed);
     alljet3eta->GetXaxis()->SetTitle("#eta");
     alljet3eta->GetYaxis()->SetTitle("Events");   
     badjet3eta->SetStats(0);
     badjet3eta->Draw("same");
     badjet3eta->SetLineColor(kBlue+1);
     TH1D *jet3etaRatio = (TH1D*)badjet3eta->Clone("jet3etaRatio");
     jet3etaRatio->SetLineColor(kBlack);
     jet3etaRatio->Sumw2();
     jet3etaRatio->SetStats(0);
     jet3etaRatio->Divide(alljet3eta);
     jet3etaRatio->Scale(100);
     goodjet3eta->Draw("same");
     goodjet3eta->SetLineColor(kGreen);
     mycanvas->SaveAs("AKTjet3eta.png");
     jet3etaRatio->Draw("HIST");
     mycanvas->SaveAs("jet3etaRatio.png");
 
     alljet4eta->Draw();
     alljet4eta->SetLineColor(kRed);
     alljet4eta->GetXaxis()->SetTitle("#eta");
     alljet4eta->GetYaxis()->SetTitle("Events");   
     badjet4eta->SetStats(0);
     badjet4eta->Draw("same");
     badjet4eta->SetLineColor(kBlue+1);
     TH1D *jet4etaRatio = (TH1D*)badjet4eta->Clone("jet4etaRatio");
     jet4etaRatio->SetLineColor(kBlack);
     jet4etaRatio->Sumw2();
     jet4etaRatio->SetStats(0);
     jet4etaRatio->Divide(alljet4eta);
     jet4etaRatio->Scale(100);
     goodjet4eta->Draw("same");
     goodjet4eta->SetLineColor(kGreen);
     mycanvas->SaveAs("AKTjet4eta.png");
     jet4etaRatio->Draw("HIST");
     mycanvas->SaveAs("jet4etaRatio.png");
*/
/*
     badjet1eta->GetXaxis()->SetTitle("#eta");
     badjet1eta->GetYaxis()->SetTitle("Events");
     badjet1eta->Draw();//"COLZ");
     mycanvas->SaveAs("badjet1eta.png");
 
     badjet2eta->GetXaxis()->SetTitle("#eta");
     badjet2eta->GetYaxis()->SetTitle("Events");
     badjet2eta->Draw();//"COLZ");
     mycanvas->SaveAs("badjet2eta.png");
 
     badjet3eta->GetXaxis()->SetTitle("#eta");
     badjet3eta->GetYaxis()->SetTitle("Events");
     badjet3eta->Draw();//"COLZ");
     mycanvas->SaveAs("badjet3eta.png");
 
     badjet4eta->GetXaxis()->SetTitle("#eta");
     badjet4eta->GetYaxis()->SetTitle("Events");
     badjet4eta->Draw();//"COLZ");
     mycanvas->SaveAs("badjet4eta.png");
*/
/* 
     badjet1theta->GetXaxis()->SetTitle("#theta");
     badjet1theta->GetYaxis()->SetTitle("Events");
     badjet1theta->Draw();//"COLZ");
     mycanvas->SaveAs("badjet1theta.png");
  
     badjet2theta->GetXaxis()->SetTitle("#theta");
     badjet2theta->GetYaxis()->SetTitle("Events");
     badjet2theta->Draw();//"COLZ");
     mycanvas->SaveAs("badjet2theta.png");
  
     badjet3theta->GetXaxis()->SetTitle("#theta");
     badjet3theta->GetYaxis()->SetTitle("Events");
     badjet3theta->Draw();//"COLZ");
     mycanvas->SaveAs("badjet3theta.png");
  
     badjet4theta->GetXaxis()->SetTitle("#theta");
     badjet4theta->GetYaxis()->SetTitle("Events");
     badjet4theta->Draw();//"COLZ");
     mycanvas->SaveAs("badjet4theta.png");
*/
     alljet1theta->Draw();
     alljet1theta->SetLineColor(kRed);
     alljet1theta->GetXaxis()->SetTitle("#theta");
     alljet1theta->GetYaxis()->SetTitle("Events");   
     badjet1theta->SetStats(0);
     badjet1theta->Draw("same");
     badjet1theta->SetLineColor(kBlue+1);
     TH1D *jet1thetaRatio = (TH1D*)badjet1theta->Clone("jet1thetaRatio");
     jet1thetaRatio->SetLineColor(kBlack);
     jet1thetaRatio->Sumw2();
     jet1thetaRatio->SetStats(0);
     jet1thetaRatio->Divide(alljet1theta);
     jet1thetaRatio->Scale(100);
     goodjet1theta->Draw("same");
     goodjet1theta->SetLineColor(kGreen);
     mycanvas->SaveAs("AKTjet1theta.png");
     jet1thetaRatio->Draw("HIST");
     mycanvas->SaveAs("jet1thetaRatio.png");
 
     alljet2theta->Draw();
     alljet2theta->SetLineColor(kRed);
     alljet2theta->GetXaxis()->SetTitle("#theta");
     alljet2theta->GetYaxis()->SetTitle("Events");   
     badjet2theta->SetStats(0);
     badjet2theta->Draw("same");
     badjet2theta->SetLineColor(kBlue+1);
     TH1D *jet2thetaRatio = (TH1D*)badjet2theta->Clone("jet2thetaRatio");
     jet2thetaRatio->SetLineColor(kBlack);
     jet2thetaRatio->Sumw2();
     jet2thetaRatio->SetStats(0);
     jet2thetaRatio->Divide(alljet2theta);
     jet2thetaRatio->Scale(100);
     goodjet2theta->Draw("same");
     goodjet2theta->SetLineColor(kGreen);
     mycanvas->SaveAs("AKTjet2theta.png");
     jet2thetaRatio->Draw("HIST");
     mycanvas->SaveAs("jet2thetaRatio.png");
 
     alljet3theta->Draw();
     alljet3theta->SetLineColor(kRed);
     alljet3theta->GetXaxis()->SetTitle("#theta");
     alljet3theta->GetYaxis()->SetTitle("Events");   
     badjet3theta->SetStats(0);
     badjet3theta->Draw("same");
     badjet3theta->SetLineColor(kBlue+1);
     TH1D *jet3thetaRatio = (TH1D*)badjet3theta->Clone("jet3thetaRatio");
     jet3thetaRatio->SetLineColor(kBlack);
     jet3thetaRatio->Sumw2();
     jet3thetaRatio->SetStats(0);
     jet3thetaRatio->Divide(alljet3theta);
     jet3thetaRatio->Scale(100);
     goodjet3theta->Draw("same");
     goodjet3theta->SetLineColor(kGreen);
     mycanvas->SaveAs("AKTjet3theta.png");
     jet3thetaRatio->Draw("HIST");
     mycanvas->SaveAs("jet3thetaRatio.png");
 

     alljet4theta->Draw();
     alljet4theta->SetLineColor(kRed);
     alljet4theta->GetXaxis()->SetTitle("#theta");
     alljet4theta->GetYaxis()->SetTitle("Events");   
     badjet4theta->SetStats(0);
     badjet4theta->Draw("same");
     badjet4theta->SetLineColor(kBlue+1);
     TH1D *jet4thetaRatio = (TH1D*)badjet4theta->Clone("jet4thetaRatio");
     jet4thetaRatio->SetLineColor(kBlack);
     jet4thetaRatio->Sumw2();
     jet4thetaRatio->SetStats(0);
     jet4thetaRatio->Divide(alljet4theta);
     jet4thetaRatio->Scale(100);
     goodjet4theta->Draw("same");
     goodjet4theta->SetLineColor(kGreen);
     mycanvas->SaveAs("AKTjet4theta.png");
     jet4thetaRatio->Draw("HIST");
     mycanvas->SaveAs("jet4thetaRatio.png");
       
// Plot DeltaR distribution
     alljet1DeltaR->Draw();
     alljet1DeltaR->SetLineColor(kRed);
     badjet1DeltaR->GetXaxis()->SetTitle("#Delta_{R}");
     badjet1DeltaR->GetYaxis()->SetTitle("Events");   
     badjet1DeltaR->SetStats(0);
     badjet1DeltaR->Draw("same");
     badjet1DeltaR->SetLineColor(kBlue+1);
     TH1D *jet1DeltaRratio = (TH1D*)badjet1DeltaR->Clone("jet1DeltaRratio");
     jet1DeltaRratio->SetLineColor(kBlack);
     jet1DeltaRratio->Sumw2();
     jet1DeltaRratio->SetStats(0);
     jet1DeltaRratio->Divide(alljet1DeltaR);
     jet1DeltaRratio->Scale(100);
     jet1DeltaRratio->Draw("sameHIST");
     mycanvas->SaveAs("badjet1DeltaR.png");
 
     alljet2DeltaR->Draw();
     alljet2DeltaR->SetLineColor(kRed);
     badjet2DeltaR->GetXaxis()->SetTitle("#Delta_{R}");
     badjet2DeltaR->GetYaxis()->SetTitle("Events");   
     badjet2DeltaR->SetStats(0);
     badjet2DeltaR->Draw("same");
     badjet2DeltaR->SetLineColor(kBlue+1);
     TH1D *jet2DeltaRratio = (TH1D*)badjet2DeltaR->Clone("jet2DeltaRratio");
     jet2DeltaRratio->SetLineColor(kBlack);
     jet2DeltaRratio->Sumw2();
     jet2DeltaRratio->SetStats(0);
     jet2DeltaRratio->Divide(alljet2DeltaR);
     jet2DeltaRratio->Scale(100);
     jet2DeltaRratio->Draw("sameHIST");
     mycanvas->SaveAs("badjet2DeltaR.png");
 
     alljet3DeltaR->Draw();
     alljet3DeltaR->SetLineColor(kRed);
     badjet3DeltaR->GetXaxis()->SetTitle("#Delta_{R}");
     badjet3DeltaR->GetYaxis()->SetTitle("Events");   
     badjet3DeltaR->SetStats(0);
     badjet3DeltaR->Draw("same");
     badjet3DeltaR->SetLineColor(kBlue+1);
     TH1D *jet3DeltaRratio = (TH1D*)badjet3DeltaR->Clone("jet3DeltaRratio");
     jet3DeltaRratio->SetLineColor(kBlack);
     jet3DeltaRratio->Sumw2();
     jet3DeltaRratio->SetStats(0);
     jet3DeltaRratio->Divide(alljet3DeltaR);
     jet3DeltaRratio->Scale(100);
     jet3DeltaRratio->Draw("sameHIST");
     mycanvas->SaveAs("badjet3DeltaR.png");
 
     alljet4DeltaR->Draw();
     alljet4DeltaR->SetLineColor(kRed);
     badjet4DeltaR->GetXaxis()->SetTitle("#Delta_{R}");
     badjet4DeltaR->GetYaxis()->SetTitle("Events");   
     badjet4DeltaR->SetStats(0);
     badjet4DeltaR->Draw("same");
     badjet4DeltaR->SetLineColor(kBlue+1);
     TH1D *jet4DeltaRratio = (TH1D*)badjet4DeltaR->Clone("jet4DeltaRratio");
     jet4DeltaRratio->SetLineColor(kBlack);
     jet4DeltaRratio->Sumw2();
     jet4DeltaRratio->SetStats(0);
     jet4DeltaRratio->Divide(alljet4DeltaR);
     jet4DeltaRratio->Scale(100);
     jet4DeltaRratio->Draw("sameHIST");
     mycanvas->SaveAs("badjet4DeltaR.png");
     cout <<endl<< "Output in TTree..."<<endl;

     AKTjetMass1->Write();
     AKTjetMass2->Write();
     //GenAKTMass2->Write();
     AKTGenMass1Comp->Write();
     AKTGenMass2Comp->Write();
     //GenUncutMass2->Write();
     AKTGenPt1Comp->Write();
     AKTGenPt2Comp->Write();
     AKTGenPt3Comp->Write();
     AKTGenPt4Comp->Write();
    
     AKTjetPT_Theta->Write();
     GenJetPT_Theta->Write();

     jet1Reso_Pt->Write();
     jet2Reso_Pt->Write();
     jet3Reso_Pt->Write();
     jet4Reso_Pt->Write();

     jet1Reso_DeltaR->Write();
     jet2Reso_DeltaR->Write();
     jet3Reso_DeltaR->Write();
     jet4Reso_DeltaR->Write();

     badjet1DeltaR->Write();
     badjet2DeltaR->Write();
     badjet3DeltaR->Write();
     badjet4DeltaR->Write();
     alljet1DeltaR->Write();
     alljet2DeltaR->Write();
     alljet3DeltaR->Write();
     alljet4DeltaR->Write();
/*
     badjet1eta->Write();
     badjet2eta->Write();
     badjet3eta->Write();
     badjet4eta->Write();
*/
     badjet1theta->Write();
     badjet2theta->Write();
     badjet3theta->Write();
     badjet4theta->Write();

     goodjet1theta->Write();
     goodjet2theta->Write();
     goodjet3theta->Write();
     goodjet4theta->Write();

     alljet1theta->Write();
     alljet2theta->Write();
     alljet3theta->Write();
     alljet4theta->Write();

     tree_output->Write();

     output->Close();
     file_sig->Close();
}

