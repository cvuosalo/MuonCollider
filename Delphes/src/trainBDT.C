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

void trainBDT(const char *BDToutputFileName) {
    TMVA::Tools::Instance();
    gROOT -> ProcessLine("gErrorIgnoreLevel = kWarning;");
    gErrorIgnoreLevel = kWarning;
    TFile* BDToutputFile = TFile::Open(BDToutputFileName, "RECREATE");
    TMVA::Factory *factory = new TMVA::Factory("TMVAClassification",BDToutputFile,"V:!Silent:Color:Transformations=I:DrawProgressBar:AnalysisType=Classification"); 
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
/* 
    dataloader -> AddVariable("BDTjet1pt1",'F');
    dataloader -> AddVariable("BDTjet1pt2",'F');
    dataloader -> AddVariable("BDTjet2pt1",'F');
    dataloader -> AddVariable("BDTjet2pt2",'F');
*/ 
    dataloader -> AddVariable("BDTjet1eta1",'F');
    dataloader -> AddVariable("BDTjet1eta2",'F');
    dataloader -> AddVariable("BDTjet2eta1",'F');
    dataloader -> AddVariable("BDTjet2eta2",'F');
/*   
    dataloader -> AddVariable("BDTjet1phi1",'F');
    dataloader -> AddVariable("BDTjet1phi2",'F');
    dataloader -> AddVariable("BDTjet2phi1",'F');
    dataloader -> AddVariable("BDTjet2phi2",'F');
  
    dataloader -> AddVariable("BDTjet1BTag1",'I');
    dataloader -> AddVariable("BDTjet1BTag2",'I');
    dataloader -> AddVariable("BDTjet2BTag1",'I');
    dataloader -> AddVariable("BDTjet2BTag2",'I');
*/     
    //dataloader -> AddVariable("BDThiggs1pt",'F');
    //dataloader -> AddVariable("BDThiggs1eta",'F');
    //dataloader -> AddVariable("BDThiggs1phi",'F');
    dataloader -> AddVariable("BDThiggs1invm",'F');

    //dataloader -> AddVariable("BDThiggs2pt",'F');
    //dataloader -> AddVariable("BDThiggs2eta",'F');
    //dataloader -> AddVariable("BDThiggs2phi",'F');
    dataloader -> AddVariable("BDThiggs2invm",'F');

    dataloader -> AddVariable("BDTdihiggspt",'F');
    //dataloader -> AddVariable("BDTdihiggseta",'F');
    //dataloader -> AddVariable("BDTdihiggsphi",'F');
    //dataloader -> AddVariable("BDTdihiggsinvm",'F');

    TCut mycuts = TCut("");
    TCut mycutb = TCut("");

    dataloader -> PrepareTrainingAndTestTree(mycuts,mycutb,"random"); 
    factory -> BookMethod(dataloader, TMVA::Types::kBDT,"BDTG","!H:!V:NTrees=3000:MinNodeSize=1.5%:BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:BaggedSampleFraction=0.6:nCuts=20:MaxDepth=4"); 
    factory -> BookMethod(dataloader, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=sigmoid:VarTransform=N:NCycles=5000:HiddenLayers=16,4:TestRate=10:LearningRate=0.02:!UseRegulator" );
    factory -> BookMethod(dataloader, TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=16,4:LearningMethod=BFGS:ValidationFraction=0.3");

    //factory -> BookMethod(dataloader, TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10");

    factory -> TrainAllMethods();
    factory -> TestAllMethods();
    factory -> EvaluateAllMethods();
    BDToutputFile  ->  Close();
    dataFile  ->  Close();
    cout << endl << "AUC for MLP method is: " << factory -> GetROCIntegral(dataloader, "MLP") << endl;
    TMVA::TMVAGui(BDToutputFileName);
    gApplication  ->  Run();
}
