void singleOutFile_maker() {
    
  // List of input files for each tree
  vector<vector<TString>> inputFiles = {
   {"/lustrehome/azaza/HccAnalysis/CMSSW_12_4_3/src/BDT_MVA/Output_files/rootFiles/output_data_2023C_Hcc.root"}, // Mass_and_BDT_DATA
   {"/lustrehome/azaza/HccAnalysis/CMSSW_12_4_3/src/BDT_MVA/Output_files/rootFiles/output_MC_Zqq2j_pt-100-200_Hcc.root", "/lustrehome/azaza/HccAnalysis/CMSSW_12_4_3/src/BDT_MVA/Output_files/rootFiles/output_MC_Zqq2j_pt-200-400_Hcc.root", "/lustrehome/azaza/HccAnalysis/CMSSW_12_4_3/src/BDT_MVA/Output_files/rootFiles/output_MC_Zqq2j_pt-400-600_Hcc.root", "/lustrehome/azaza/HccAnalysis/CMSSW_12_4_3/src/BDT_MVA/Output_files/rootFiles/output_MC_Zqq2j_pt-600_Hcc.root"}, // Mass_and_BDT_ZJets
   {"/lustrehome/azaza/HccAnalysis/CMSSW_12_4_3/src/BDT_MVA/Output_files/rootFiles/output_MC_VBFZqq_Hcc.root"}, // Mass_and_BDT_ZJets_EWK
   {"/lustrehome/azaza/HccAnalysis/CMSSW_12_4_3/src/BDT_MVA/Output_files/rootFiles/output_MC_Wqq2j_pt-100-200_Hcc.root","/lustrehome/azaza/HccAnalysis/CMSSW_12_4_3/src/BDT_MVA/Output_files/rootFiles/output_MC_Wqq2j_pt-200-400_Hcc.root","/lustrehome/azaza/HccAnalysis/CMSSW_12_4_3/src/BDT_MVA/Output_files/rootFiles/output_MC_Wqq2j_pt-400-600_Hcc.root","/lustrehome/azaza/HccAnalysis/CMSSW_12_4_3/src/BDT_MVA/Output_files/rootFiles/output_MC_Wqq2j_pt-600_Hcc.root"}, //Mass_and_BDT_WJets
   {"/lustrehome/azaza/HccAnalysis/CMSSW_12_4_3/src/BDT_MVA/Output_files/rootFiles/output_MC_VBFWqq_Hcc.root"}, //Mass_and_BDT_WJets_EWK
   {"/lustrehome/azaza/HccAnalysis/CMSSW_12_4_3/src/BDT_MVA/Output_files/rootFiles/output_MC_VBFHCC_Hcc.root"}, // Mass_and_BDT_VBF_Hcc_Dipole
   {"/lustrehome/azaza/HccAnalysis/CMSSW_12_4_3/src/BDT_MVA/Output_files/rootFiles/output_MC_ggHCC_Hcc.root"}, //Mass_and_BDT_ggF_Hcc
   {"/lustrehome/azaza/HccAnalysis/CMSSW_12_4_3/src/BDT_MVA/Output_files/rootFiles/output_MC_VBFHBB_Hcc.root"}, //Mass_and_BDT_VBF_Hbb
   {"/lustrehome/azaza/HccAnalysis/CMSSW_12_4_3/src/BDT_MVA/Output_files/rootFiles/output_MC_ggHBB_Hcc.root"}, //Mass_and_BDT_ggF_Hbb
   };
  vector<TString> treeNames = {"Mass_and_BDT_DATA", "Mass_and_BDT_ZJets", "Mass_and_BDT_ZJets_EWK", "Mass_and_BDT_WJets", "Mass_and_BDT_WJets_EWK", "Mass_and_BDT_VBF_Hcc_Dipol", "Mass_and_BDT_ggF_Hcc", "Mass_and_BDT_VBF_Hbb_Dipol", "Mass_and_BDT_ggF_Hbb" };

  TString outputFileName = "mcc_and_bdt_all_Nom.root";
  gSystem->Exec("rm -f " + outputFileName); // Remove existing output file if it exists
  TFile *outputFile = TFile::Open(outputFileName, "RECREATE");   

  const int batchSize = 1000; // Set an appropriate batch size

  double mass_cJets, BDT_score, xs, lumi;

  for (size_t treeIndex = 0; treeIndex < inputFiles.size(); ++treeIndex) {
       const vector<TString> &filesForTree = inputFiles[treeIndex];
       const TString &treeName = treeNames[treeIndex];

    TNtuple *tree = new TNtuple(treeName, "InvariantMass", "mcc:bdtout:weight");

    for (const TString &inputFile : filesForTree) {
      TFile *f_bkg = TFile::Open(inputFile);
      TH1* histo_counts= dynamic_cast<TH1*>(f_bkg->Get("CutEff_NEvents"));
      TTree *t_bkg = (TTree *)f_bkg->Get("T");
      
      t_bkg->SetBranchAddress("mass_cJets", &mass_cJets);
      t_bkg->SetBranchAddress("BDT_score", &BDT_score);
      t_bkg->SetBranchAddress("xs", &xs);
      t_bkg->SetBranchAddress("lumi", &lumi);
      

      for (Long64_t ievt = 0; ievt < t_bkg->GetEntries(); ievt++) {
        t_bkg->GetEntry(ievt);
        lumi=17.7799;
        float weight = xs*lumi/(histo_counts->GetBinContent(1));
        tree->Fill(mass_cJets,BDT_score,weight);
      }

      f_bkg->Close();
    }
    
    outputFile->cd();
    tree->Write();
    delete tree;  
  }


  outputFile->Close();
  delete outputFile;
}
