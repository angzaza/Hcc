#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/CrossValidation.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

using namespace TMVA;

void TMVA_Classification(){

	// Output file
	TFile *fout = new TFile("TMVA_output_v1.root", "RECREATE");
	
	// SIGNAL tree	
	TFile *f_sig = new TFile("AnalysedTree_MC_VBFHCC_Hcc_merged_v2.root");
	TTree *signalTree     = (TTree*)f_sig->Get("FinalTree");

	//BKG tree
	TFile *f_bkg = new TFile("AnalysedTree_data_2023C_Hcc_merged_v2.root");
	TTree *bkgTree     = (TTree*)f_bkg->Get("FinalTree");
	
	//weights
	Double_t sigWeight1  = 1.0;
  Double_t bkgWeight1 = 1.0;
    
    
  Factory *factory = new Factory("TMVA_new", fout, "");
  DataLoader *dataloader = new DataLoader("dataset_v1");
    	
  dataloader->AddSignalTree(signalTree,sigWeight1);
  dataloader->AddBackgroundTree(bkgTree,bkgWeight1);
     	

  dataloader->AddVariable("CvsB_jetC1", 'D');
  dataloader->AddVariable("CvsB_jetC2", 'D');
  dataloader->AddVariable("CvsL_jetC1", 'D');
  dataloader->AddVariable("CvsL_jetC2", 'D');
  dataloader->AddVariable("mqq", 'D');
  dataloader->AddVariable("Deta_qq", 'D');
  dataloader->AddVariable("Dphi_qq", 'D');
  dataloader->AddVariable("Alpha_qq", 'D');
  //dataloader->AddVariable("qgl_VBF1", 'D');
  //dataloader->AddVariable("qgl_VBF2", 'D');
  dataloader->AddVariable("pz_4jets", 'D');
  dataloader->AddVariable("pt_norm", 'D');
  dataloader->AddVariable("DR_HiggsVBF1", 'D');
  dataloader->AddVariable("DR_HiggsVBF2", 'D');
  dataloader->AddVariable("Dphi_qq_cc", 'D');
  dataloader->AddVariable("njets", 'I');
  dataloader->AddVariable("jetEne_sum", 'D');
  dataloader->AddVariable("jetPt_sum", 'D');
    	
    	//cuts
    	TCut cutS="";
    	TCut cutB="";
    	
    	dataloader->PrepareTrainingAndTestTree( cutS, cutB,"SplitMode=Random:NormMode=NumEvents:!V" );
	
	// Booking of MVA methods : BDT
        //factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",
        //                   "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=50" );

	factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG",
											"!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:PruneMethod=CostComplexity:PruneStrength=50:NNodesMax=5" );
        // Training the MVA methods
        factory->TrainAllMethods();
        
        // Testing the MVA methods
        factory->TestAllMethods();
        
        // Evaluating the MVA methods
        factory->EvaluateAllMethods();
		
	// Save the output
    fout->Close();
    
    std::cout << "==> Wrote root file: " << fout->GetName() << std::endl;
    std::cout << "==> TMVAClassification is done!" << std::endl;
    
    delete factory;
    delete dataloader;
    
    // Launch the GUI for the root macros
    if (!gROOT->IsBatch()){
        TMVAGui("TMVA_output_v1.root");
    	}


}
