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

#include "TMVA/Reader.h"

#include "TMVA/Tools.h"
using namespace TMVA;
#endif
void applyBDTindividual(TTree* theTree, TH1F* histBdt, TH1F* jetpair1, TH1F* jetpair2){
    TMVA::Tools::Instance();
    //gROOT -> ProcessLine("gErrorIgnoreLevel = kWarning;");
    //gErrorIgnoreLevel = kWarning;
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

    Float_t BDTjet1BTag1;
    Float_t BDTjet1BTag2;
    Float_t BDTjet2BTag1;
    Float_t BDTjet2BTag2;

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

    Float_t NBTag;
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

    Int_t BDTjet1BTag1_d;
    Int_t BDTjet1BTag2_d;
    Int_t BDTjet2BTag1_d;
    Int_t BDTjet2BTag2_d;

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
   
    reader -> AddVariable("BDTjet1BTag1", &BDTjet1BTag1);
    reader -> AddVariable("BDTjet1BTag2", &BDTjet1BTag2);
    reader -> AddVariable("BDTjet2BTag1", &BDTjet2BTag1);
    reader -> AddVariable("BDTjet2BTag2", &BDTjet2BTag2);

    reader -> AddVariable("BDTMissingETMET", &BDTMissingETMET);
    reader -> AddVariable("BDTMissingETEta", &BDTMissingETEta);
    reader -> AddVariable("BDTMissingETPhi", &BDTMissingETPhi);
  
    reader -> AddVariable("NBTag", &NBTag);

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
    /*
    theTree -> SetBranchAddress("BDTNjets", &BDTNjets_d); 
    
    theTree -> SetBranchAddress("BDTjet1pt1", &BDTjet1pt1_d);
    theTree -> SetBranchAddress("BDTjet1pt2", &BDTjet1pt2_d);
    theTree -> SetBranchAddress("BDTjet2pt1", &BDTjet2pt1_d);
    theTree -> SetBranchAddress("BDTjet2pt2", &BDTjet2pt2_d);

    theTree -> SetBranchAddress("BDTjet1eta1", &BDTjet1eta1_d);
    theTree -> SetBranchAddress("BDTjet1eta2", &BDTjet1eta2_d);
    theTree -> SetBranchAddress("BDTjet2eta1", &BDTjet2eta1_d);
    theTree -> SetBranchAddress("BDTjet2eta2", &BDTjet2eta2_d);

    theTree -> SetBranchAddress("BDTjet1phi1", &BDTjet1phi1_d);
    theTree -> SetBranchAddress("BDTjet1phi2", &BDTjet1phi2_d);
    theTree -> SetBranchAddress("BDTjet2phi1", &BDTjet2phi1_d);
    theTree -> SetBranchAddress("BDTjet2phi2", &BDTjet2phi2_d);

    theTree -> SetBranchAddress("BDTjet1BTag1", &BDTjet1BTag1_d);
    theTree -> SetBranchAddress("BDTjet1BTag2", &BDTjet1BTag2_d);
    theTree -> SetBranchAddress("BDTjet2BTag1", &BDTjet2BTag1_d);
    theTree -> SetBranchAddress("BDTjet2BTag2", &BDTjet2BTag2_d);
    
    theTree -> SetBranchAddress("NBTag", &NBTag_d);

    theTree -> SetBranchAddress("BDTMissingETMET", &BDTMissingETMET_d);
    theTree -> SetBranchAddress("BDTMissingETEta", &BDTMissingETEta_d);
    theTree -> SetBranchAddress("BDTMissingETPhi", &BDTMissingETPhi_d);
   
    theTree -> SetBranchAddress("DiHiggsMETDeltaTheta", &DiHiggsMETDeltaTheta_d);

    theTree -> SetBranchAddress("BDThiggsdeltaPhi", &BDThiggsdeltaPhi_d);
    theTree -> SetBranchAddress("BDThiggsdeltaEta", &BDThiggsdeltaEta_d);

    theTree -> SetBranchAddress("BDThiggs1pt", &BDThiggs1pt_d);
    theTree -> SetBranchAddress("BDThiggs1eta", &BDThiggs1eta_d);
    theTree -> SetBranchAddress("BDThiggs1invm", &BDThiggs1invm_d);

    theTree -> SetBranchAddress("BDThiggs2pt", &BDThiggs2pt_d);
    theTree -> SetBranchAddress("BDThiggs2eta", &BDThiggs2eta_d);
    theTree -> SetBranchAddress("BDThiggs2invm", &BDThiggs2invm_d);

    theTree -> SetBranchAddress("BDTdihiggspt", &BDTdihiggspt_d);
    theTree -> SetBranchAddress("BDTdihiggseta", &BDTdihiggseta_d);
    theTree -> SetBranchAddress("BDTdihiggsinvm", &BDTdihiggsinvm_d);
    */

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

    TLeaf *jet1BTag1 = theTree -> GetLeaf("BDTjet1BTag1");
    TLeaf *jet1BTag2 = theTree -> GetLeaf("BDTjet1BTag2");
    TLeaf *jet2BTag1 = theTree -> GetLeaf("BDTjet2BTag1");
    TLeaf *jet2BTag2 = theTree -> GetLeaf("BDTjet2BTag2");
    
    TLeaf *NBTagLeaf = theTree -> GetLeaf("NBTag");
    
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
/*
	BDTNjets = BDTNjets_d;

        BDTjet1pt1 = (Float_t) BDTjet1pt1_d;
	BDTjet1pt2 = (Float_t) BDTjet1pt2_d;
	BDTjet2pt1 = (Float_t) BDTjet2pt1_d;
	BDTjet2pt2 = (Float_t) BDTjet2pt2_d;
	
	BDTjet1eta1 = (Float_t) BDTjet1eta1_d;
	BDTjet1eta2 = (Float_t) BDTjet1eta2_d;
	BDTjet2eta1 = (Float_t) BDTjet2eta1_d;
	BDTjet2eta2 = (Float_t) BDTjet2eta2_d;
	
	BDTjet1phi1 = (Float_t) BDTjet1phi1_d;
	BDTjet1phi2 = (Float_t) BDTjet1phi2_d;
	BDTjet2phi1 = (Float_t) BDTjet2phi1_d;
	BDTjet2phi2 = (Float_t) BDTjet2phi2_d;

	BDTjet1BTag1 = (Float_t) BDTjet1BTag1_d;
	BDTjet1BTag2 = (Float_t) BDTjet1BTag2_d;
	BDTjet2BTag1 = (Float_t) BDTjet2BTag1_d;
	BDTjet2BTag2 = (Float_t) BDTjet2BTag2_d;
	
	BDTMissingETMET = (Float_t) BDTMissingETMET_d;
	BDTMissingETEta = (Float_t) BDTMissingETEta_d;
	BDTMissingETPhi = (Float_t) BDTMissingETPhi_d;

	BDThiggsdeltaEta = 0;
	BDThiggsdeltaPhi = (Float_t) BDThiggsdeltaPhi_d;


	BDThiggs1pt = (Float_t) BDThiggs1pt_d;
	BDThiggs1eta = (Float_t) BDThiggs1eta_d;
	BDThiggs1invm = (Float_t) BDThiggs1invm_d;
	
	BDThiggs2pt = (Float_t) BDThiggs2pt_d;
	BDThiggs2eta = (Float_t) BDThiggs2eta_d;
	BDThiggs2invm = (Float_t) BDThiggs2invm_d;
	
	BDTdihiggspt = (Float_t) BDTdihiggspt_d;
	BDTdihiggseta = (Float_t) BDTdihiggseta_d;
	BDTdihiggsinvm = (Float_t) BDTdihiggsinvm_d;
	
	NBTag = (Float_t) NBTag_d;
	DiHiggsMETDeltaTheta = (Float_t) DiHiggsMETDeltaTheta_d;
*/	
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

	BDTjet1BTag1 = jet1BTag1 -> GetValue();
	BDTjet1BTag2 = jet1BTag2 -> GetValue();
	BDTjet2BTag1 = jet2BTag1 -> GetValue();
	BDTjet2BTag2 = jet2BTag2 -> GetValue();
	
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
	
	NBTag =  NBTagLeaf -> GetValue();
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

    THStack* JetPair1BDT = new THStack("JetPair1BDT", "");
    THStack* JetPair2BDT = new THStack("JetPair2BDT", "");
    TH1F *AKTjetMass1_sigBDT = new TH1F("AKTjetMass1_sigBDT", "Anti_KTjet leading jets pair invariant mass", 50, 60, 190);
    TH1F *AKTjetMass2_sigBDT = new TH1F("AKTjetMass2_sigBDT", "Anti_KTjet sub-leading jets pair invariant mass", 50, 60, 190); 

    TH1F *AKTjetMass1_bkg1BDT = new TH1F("AKTjetMass1_bkg1BDT", "Anti_KTjet leading jets pair invariant mass", 50, 60, 190);
    TH1F *AKTjetMass2_bkg1BDT = new TH1F("AKTjetMass2_bkg1BDT", "Anti_KTjet sub-leading jets pair invariant mass", 50, 60, 190); 

    TH1F *AKTjetMass1_bkg2BDT = new TH1F("AKTjetMass1_bkg2BDT", "Anti_KTjet leading jets pair invariant mass", 50, 60, 190);
    TH1F *AKTjetMass2_bkg2BDT = new TH1F("AKTjetMass2_bkg2BDT", "Anti_KTjet sub-leading jets pair invariant mass", 50, 60, 190); 

    TH1F *AKTjetMass1_bkg3BDT = new TH1F("AKTjetMass1_bkg3BDT", "Anti_KTjet leading jets pair invariant mass", 50, 60, 190);
    TH1F *AKTjetMass2_bkg3BDT = new TH1F("AKTjetMass2_bkg3BDT", "Anti_KTjet sub-leading jets pair invariant mass", 50, 60, 190);

    Double_t weight1 = 10*0.0008201/(0.0008201 + 0.0009227 + 0.003743 + 0.03168);
    Double_t weight2 = 10*0.0009227/(0.0008201 + 0.0009227 + 0.003743 + 0.03168);
    Double_t weight3 = 10*0.003743/(0.0008201 + 0.0009227 + 0.003743 + 0.03168);
    Double_t weight4 = 10*0.03168/(0.0008201 + 0.0009227 + 0.003743 + 0.03168);
 
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

    TCanvas *totalcanvas = new TCanvas("totalcanvas", "Canvas", 300, 300, 1200, 1200);
    gStyle->SetHistMinimumZero();
    totalcanvas -> SetLogy();
    totalcanvas -> SetWindowSize(1500,1500);
    totalcanvas -> SetCanvasSize(1200,1200);
    JetPair1BDT -> Draw("HIST");
    AKTjetMass1_sigBDT -> Draw("sameHIST");
    TLegend *legend = new TLegend(0.65, 0.65, 0.875, 0.85);
    legend -> AddEntry(AKTjetMass1_sigBDT, "signal HH", "l");
    legend -> AddEntry(AKTjetMass1_bkg1BDT, "QCD", "f");
    legend -> AddEntry(AKTjetMass1_bkg2BDT, "single Higgs", "f");
    legend -> AddEntry(AKTjetMass1_bkg3BDT, "Z + jets", "f");
    legend -> SetBorderSize(0);
    legend -> Draw();
    JetPair1BDT -> GetXaxis() -> SetTitle("m_{H_1} [GeV]");
    JetPair1BDT -> SetMaximum(1e5);
    JetPair1BDT -> GetYaxis() -> SetTitle("events");
    TLatex title1(85, 2.4e5, "#bf{#font[72]{Muon Collider Simulation (Delphes)}}");
    title1.SetTextSize(0.04); 
    TLatex latexUp1(70, 5e4, "\\sqrt{s} = 3 TeV}");
    TLatex latexDown1(75, 2e4, "#font[52]{L} #font[42]{=} #font[52]{1 ab^{-1}}");
    latexUp1.SetTextSize(0.04); 
    latexDown1.SetTextSize(0.04); 
    latexUp1.Draw("same");
    latexDown1.Draw("same");
    title1.Draw("same");
    totalcanvas -> SaveAs("JetPairs1_sig+bkgBDT.png");

    JetPair2BDT -> Draw("HIST");
    AKTjetMass2_sigBDT -> Draw("sameHIST");
    TLegend *legend2 = new TLegend(0.65, 0.65, 0.875, 0.85);
    legend2 -> AddEntry(AKTjetMass2_sigBDT, "signal HH", "l");
    legend2 -> AddEntry(AKTjetMass2_bkg1BDT, "QCD", "f");
    legend2 -> AddEntry(AKTjetMass2_bkg2BDT, "single Higgs ", "f");
    legend2 -> AddEntry(AKTjetMass2_bkg3BDT, "Z + jets", "f");
    legend2 -> SetBorderSize(0);
    legend2 -> Draw();
    JetPair2BDT -> GetXaxis() -> SetTitle("m_{H_2} [GeV]");
    JetPair2BDT -> GetYaxis() -> SetTitle("events");
    JetPair2BDT -> SetMaximum(1e5);
    TLatex title2(85, 2.5e5, "#bf{#font[72]{Muon Collider Simulation (Delphes)}}");
    title2.SetTextSize(0.04); 
    TLatex latexUp2(70, 5e4, "\\sqrt{s} = 3 TeV}");
    TLatex latexDown2(75, 1.75e4, "#font[52]{L} #font[42]{=} #font[52]{1 ab^{-1}}");
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
   
    BDTs_b-> Draw("HIST");
    histBdt_sig -> Draw("sameHIST");
    TLegend *legendBdt = new TLegend(0.65, 0.65, 0.875, 0.85);
    legendBdt -> AddEntry(histBdt_sig, "signal HH", "l");
    legendBdt -> AddEntry(histBdt_bkg1, "QCD", "f");
    legendBdt -> AddEntry(histBdt_bkg2, "single Higgs", "f");
    legendBdt -> AddEntry(histBdt_bkg3, "Z + jets", "f");
    legendBdt -> SetBorderSize(0);
    legendBdt -> Draw();
    BDTs_b -> GetXaxis() -> SetTitle("BDT_{Z+jet} score");
    BDTs_b -> GetYaxis() -> SetTitle("events");
    TLatex titleBdt(0.41, 3800, "#bf{#font[72]{Muon Collider Simulation (Delphes)}}");
    titleBdt.SetTextSize(0.04); 
    TLatex latexUpBdt(0.365, 1500, "\\sqrt{s} = 3 TeV}");
    TLatex latexDownBdt(0.375, 900, "#font[52]{L = 1 ab^{-1}}");
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

