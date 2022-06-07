//usage: root -l macroExample.C\(\"inputfile.root\"\,\"outputfile.root\"\)
#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
//put header files you need here
#endif

void macroExample(const char *inputFile, const char *outputFile){
     gSystem -> Load("libDelphes.so");
     TFile *file_sig = new TFile(inputFile);
     TFile *output = new TFile(outputFile, "recreate");
     TTree *tree_sig = (TTree*)file_sig -> Get("Delphes");
     TTree *tree_output = new TTree("tree_output","Delphes");

     Int_t nEntries = tree_sig -> GetEntries();
     //TLeaf * AKTjet_size = tree_sig -> GetLeaf("AKTjet_size");
     //TLeaf * AKTjet_pt = tree_sig -> GetLeaf("AKTjet.PT");
     for (Long64_t entry = 0; entry < nEntries; entry++){
	 tree_sig -> GetEntry(entry);
	 tree_output -> GetEntry(entry);
	 //write macro algorithm here
	 /*
	 AKTjet_size -> GetBranch() -> GetEntry(entry);
	 Int_t nAKTjet = AKTjet_size -> GetValue();
	 for (Int_t aktentry = 0; aktentry < nAKTjet; aktentry++){
	     AKTjetPt = AKTjet_pt -> GetValue(aktentry);
	 }
	 */
     }
     tree_output -> Write();
     output -> Close();
     file_sig -> Close();
}

