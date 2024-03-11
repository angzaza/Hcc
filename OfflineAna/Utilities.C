#define myAnalyzer_cxx
#define NCUTS 10

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TMVA/Reader.h>


void myAnalizer::TreeFin_Init(TTree *&tree, Double_t &isMC, Double_t &lumi, Double_t &run, Double_t &evt,int &entry, Double_t &XS, Double_t &XS_err,Double_t &Lumi, Double_t &puFactor, Double_t &pt_jetC1, Double_t &pt_jetC2, Double_t &pt_jetVBF1, Double_t &pt_jetVBF2, Double_t &eta_jetC1, Double_t &eta_jetC2, Double_t &eta_jetVBF1, Double_t &eta_jetVBF2, Double_t &CvsAll_jetC1, Double_t &CvsAll_jetC2, Double_t &CvsB_jetC1, Double_t &CvsB_jetC2, Double_t &CvsL_jetC1, Double_t &CvsL_jetC2, Double_t &mqq, Double_t &Deta_qq, Double_t &Dphi_qq, Double_t  &Alpha_qq, Double_t &qgl_VBF1, Double_t &qgl_VBF2, Double_t &QvsG_VBF1, Double_t &QvsG_VBF2, Double_t &pz_4jets, Double_t &pt_norm, Double_t &DR_HiggsVBF1, Double_t &DR_HiggsVBF2, Double_t &Dphi_qq_cc, int &njets, Double_t &jetEne_sum, Double_t  &jetPt_sum, Double_t &mCC){
		// Set tree branches
		tree->Branch("isMC", &isMC);
		tree->Branch("lumi", &lumi);
		tree->Branch("run", &run);
		tree->Branch("evt", &evt);
		tree->Branch("entry", &entry);
		tree->Branch("XS", &XS);
		tree->Branch("XS_err", &XS_err);
		tree->Branch("Lumi", &Lumi);
		tree->Branch("puFactor", &puFactor);
    tree->Branch("pt_jetC1", &pt_jetC1);
    tree->Branch("pt_jetC2", &pt_jetC2);
    tree->Branch("pt_jetVBF1", &pt_jetVBF1);
    tree->Branch("pt_jetVBF2", &pt_jetVBF2);
    tree->Branch("eta_jetC1", &eta_jetC1);
    tree->Branch("eta_jetC2", &eta_jetC2);
    tree->Branch("eta_jetVBF1", &eta_jetVBF1);
    tree->Branch("eta_jetVBF2", &eta_jetVBF2);
    tree->Branch("CvsAll_jetC1", &CvsAll_jetC1);
    tree->Branch("CvsAll_jetC2", &CvsAll_jetC2);
    tree->Branch("CvsB_jetC1", &CvsB_jetC1);
    tree->Branch("CvsB_jetC2", &CvsB_jetC2);
    tree->Branch("CvsL_jetC1", &CvsL_jetC1);
    tree->Branch("CvsL_jetC2", &CvsL_jetC2);
    tree->Branch("mqq", &mqq);
    tree->Branch("Deta_qq", &Deta_qq);
    tree->Branch("Dphi_qq", &Dphi_qq);
    tree->Branch("Alpha_qq", &Alpha_qq);
    tree->Branch("qgl_VBF1", &qgl_VBF1);
    tree->Branch("qgl_VBF2", &qgl_VBF2);
    tree->Branch("QvsG_VBF1", &QvsG_VBF1);
    tree->Branch("QvsG_VBF2", &QvsG_VBF2);
    tree->Branch("pz_4jets", &pz_4jets);
    tree->Branch("pt_norm", &pt_norm);
    tree->Branch("DR_HiggsVBF1", &DR_HiggsVBF1);
    tree->Branch("DR_HiggsVBF2", &DR_HiggsVBF2);
    tree->Branch("Dphi_qq_cc", &Dphi_qq_cc);
    tree->Branch("njets", &njets);
    tree->Branch("jetEne_sum", &jetEne_sum);
    tree->Branch("jetPt_sum", &jetPt_sum);
    tree->Branch("mCC", &mCC);
}


void myAnalizer::Fill_CutName(TString listCut[NCUTS]){
    // Init a vector of strings w/ the names of the cuts
    listCut[0] = "BeforeCuts";
    listCut[1] = "HLT_fired";
    listCut[2] = "lep_MET_veto";
    listCut[3] = "pT_cut";
    listCut[4] = "jet_vetomaps";
    listCut[5] = "hltmatch";
    listCut[6] = "2_centralj";
    listCut[7] = "2c2VBF";
    listCut[8] = "CTagCut";
    listCut[9] = "VBFCut";
  
}

void myAnalizer::Draw_CutEffCanvas(TCanvas *canv, TH1I *hist, Int_t cut[NCUTS], TString listCut[NCUTS]){
    // This function writes on the canvas the histo of the cuts efficiency
    for(int k=0; k<NCUTS; k++){
      hist->Fill(k+1, cut[k]);
      hist->GetXaxis()->SetBinLabel(k+1, listCut[k]);
    }
    //    canv->SetLogy();
    hist->DrawCopy("HIST TEXT0");
    hist->Write();
    //    canv->Write();
    //    canv->Close();
}
    //

std::pair<float, float> myAnalizer::CrossSection(TString dataset){
   float xs=0., xs_err=0.;

   // cross sections in fb////////////////
   //
   // VBFHcc
   float XS_VBFHcc = 120.8438;
   float XSerr_VBFHcc = 0.;

   // VBFHbb
   float XS_VBFHbb = 2434.432;
   float XSerr_VBFHbb = 0.;

   // ggHcc
   float XS_ggHcc = 873.082;
   float XSerr_ggHcc = 0.;

   // ggHcc
   float XS_ggHbb = 17588.48;
   float XSerr_ggHbb = 0.;

   //QCD HT binned
   float XS_QCD_HT100_200 = 2.53E+10;
   float XSerr_QCD_HT100_200 = 0.;

   float XS_QCD_HT200_400 = 1.96E+09;
   float XSerr_QCD_HT200_400 = 0.;

   float XS_QCD_HT400_600 = 9.68E+07;
   float XSerr_QCD_HT400_600 = 0.;

   float XS_QCD_HT600_800 = 1.37E+07;
   float XSerr_QCD_HT600_800 = 0.;

   float XS_QCD_HT800_1000 = 3.03E+06;
   float XSerr_QCD_HT800_1000 = 0.;

   float XS_QCD_HT1000_1200 = 8.81E+05;
   float XSerr_QCD_HT1000_1200 = 0.;

   float XS_QCD_HT1200_1500 = 3.87E+05;
   float XSerr_QCD_HT1200_1500 = 0.;

   float XS_QCD_HT1500_2000 = 1.27E+05;
   float XSerr_QCD_HT1500_2000 = 0.;

   float XS_QCD_HT2000 = 2.66E+04;
   float XSerr_QCD_HT2000 = 0.;

  // TTBar 
   float XS_TTto2L2Nu = 7.62e+05;
   float XSerr_TTto2L2Nu = 1.35e+02;

   // QCD Zqq
   float XS_Zto2Q_HT200_400 = 1.09e+06;
   float XSerr_Zto2Q_HT200_400 = 6.01e+03;

   float XS_Zto2Q_HT400_600 = 1.24e+05;
   float XSerr_Zto2Q_HT400_600 = 5.03e+02;

   float XS_Zto2Q_HT600_800 = 2.78e+04;
   float XSerr_Zto2Q_HT600_800 = 1.57e+02;

   float XS_Zto2Q_HT800 = 1.45e+04;
   float XSerr_Zto2Q_HT800 = 7.30e+01;

   // QCD Wqq
   float XS_Wto2Q_HT200_400 = 2.74e+06 ;
   float XSerr_Wto2Q_HT200_400 = 1.17e+04;

   float XS_Wto2Q_HT400_600 = 3.01e+05;
   float XSerr_Wto2Q_HT400_600 = 1.32e+03;

   float XS_Wto2Q_HT600_800 = 6.50e+04;
   float XSerr_Wto2Q_HT600_800 = 2.72e+02;

   float XS_Wto2Q_HT800 = 3.21e+04;
   float XSerr_Wto2Q_HT800 = 1.74e+02;


   // Zqq pt-binned
   float XS_Zto2Q_pt100_200=2.97E+05;
   float XSerr_Zto2Q_pt100_200=0.;
   float XS_Zto2Q_pt200_400=2.14E+04;
   float XSerr_Zto2Q_pt200_400=0.;
   float XS_Zto2Q_pt400_600=7.48E+02;
   float XSerr_Zto2Q_pt400_600=0.;
   float XS_Zto2Q_pt600=8.74E+01;
   float XSerr_Zto2Q_pt600=0.;

   // Zqq pt-binned
   float XS_Wto2Q_pt100_200=1.55E+06;
   float XSerr_Wto2Q_pt100_200=0.;
   float XS_Wto2Q_pt200_400=1.04E+05;
   float XSerr_Wto2Q_pt200_400=0.;
   float XS_Wto2Q_pt400_600=3.53E+03;
   float XSerr_Wto2Q_pt400_600=0.;
   float XS_Wto2Q_pt600=4.19E+02;
   float XSerr_Wto2Q_pt600=0.;

   // VBFZqq
   float XS_VBFZto2Q=1.37E+04;
   float XSerr_VBFZto2Q=0.;

   // VBFWqq
   float XS_VBFWto2Q=9.56E+04;
   float XSerr_VBFWto2Q=0.;

   if(dataset.Contains("VBFHCC")){
      xs = XS_VBFHcc;
      xs_err = XSerr_VBFHcc;
   }
   if(dataset.Contains("VBFHBB")){
      xs = XS_VBFHbb;
      xs_err = XSerr_VBFHbb;
   }
   if(dataset.Contains("ggHCC")){
      xs = XS_ggHcc;
      xs_err = XSerr_ggHcc;
   }
   if(dataset.Contains("ggHBB")){
      xs = XS_ggHbb;
      xs_err = XSerr_ggHbb;
   }
   if(dataset.Contains("QCD-4Jets_HT-100to200")){
      xs = XS_QCD_HT100_200;
      xs_err = XSerr_QCD_HT100_200;
   }
   if(dataset.Contains("QCD-4Jets_HT-200to400")){
      xs = XS_QCD_HT200_400;
      xs_err = XSerr_QCD_HT200_400;
   }
   if(dataset.Contains("QCD-4Jets_HT-400to600")){
      xs = XS_QCD_HT400_600;
      xs_err = XSerr_QCD_HT400_600;
   }
   if(dataset.Contains("QCD-4Jets_HT-600to800")){
      xs = XS_QCD_HT600_800;
      xs_err = XSerr_QCD_HT600_800;
   }
   if(dataset.Contains("QCD-4Jets_HT-800to1000")){
      xs = XS_QCD_HT800_1000;
      xs_err = XSerr_QCD_HT800_1000;
   }
   if(dataset.Contains("QCD-4Jets_HT-1000to1200")){
      xs = XS_QCD_HT1000_1200;
      xs_err = XSerr_QCD_HT1000_1200;
   }
   if(dataset.Contains("QCD-4Jets_HT-1200to1500")){
      xs = XS_QCD_HT1200_1500;
      xs_err = XSerr_QCD_HT1200_1500;
   }
   if(dataset.Contains("QCD-4Jets_HT-1500to2000")){
      xs = XS_QCD_HT1500_2000;
      xs_err = XSerr_QCD_HT1500_2000;
   }
   if(dataset.Contains("QCD-4Jets_HT-2000")){
      xs = XS_QCD_HT2000;
      xs_err = XSerr_QCD_HT2000;
   }
   if(dataset.Contains("TTto2L2Nu")){
      xs = XS_TTto2L2Nu;
      xs_err = XSerr_TTto2L2Nu;
   }
   if(dataset.Contains("Wqq_HT-200-400")){
      xs = XS_Wto2Q_HT200_400;
      xs_err = XSerr_Wto2Q_HT200_400;
   }
   if(dataset.Contains("Wqq_HT-400-600")){
      xs = XS_Wto2Q_HT400_600;
      xs_err = XSerr_Wto2Q_HT400_600;
   }
   if(dataset.Contains("Wqq_HT-600-800")){
      xs = XS_Wto2Q_HT600_800;
      xs_err = XSerr_Wto2Q_HT600_800;
   }
   if(dataset.Contains("Wqq_HT-800")){
      xs = XS_Wto2Q_HT800;
      xs_err = XSerr_Wto2Q_HT800;
   }
   if(dataset.Contains("Zqq_HT-200-400")){
      xs = XS_Zto2Q_HT200_400;
      xs_err = XSerr_Zto2Q_HT200_400;
   }
   if(dataset.Contains("Zqq_HT-400-600")){
      xs = XS_Zto2Q_HT400_600;
      xs_err = XSerr_Zto2Q_HT400_600;
   }
   if(dataset.Contains("Zqq_HT-600-800")){
      xs = XS_Zto2Q_HT600_800;
      xs_err = XSerr_Zto2Q_HT600_800;
   }
   if(dataset.Contains("Zqq_HT-800")){
      xs = XS_Zto2Q_HT800;
      xs_err = XSerr_Zto2Q_HT800;
   }
   if(dataset.Contains("Zqq_pt-100-200")){
      xs = XS_Zto2Q_pt100_200;
      xs_err = XSerr_Zto2Q_pt100_200;
   }
   if(dataset.Contains("Zqq_pt-200-400")){
      xs = XS_Zto2Q_pt200_400;
      xs_err = XSerr_Zto2Q_pt200_400;
   }
   if(dataset.Contains("Zqq_pt-400-600")){
      xs = XS_Zto2Q_pt400_600;
      xs_err = XSerr_Zto2Q_pt400_600;
   }
   if(dataset.Contains("Zqq_pt-600")){
      xs = XS_Zto2Q_pt600;
      xs_err = XSerr_Zto2Q_pt600;
   }
   if(dataset.Contains("Wqq_pt-100-200")){
      xs = XS_Wto2Q_pt100_200;
      xs_err = XSerr_Wto2Q_pt100_200;
   }
   if(dataset.Contains("Wqq_pt-200-400")){
      xs = XS_Wto2Q_pt200_400;
      xs_err = XSerr_Wto2Q_pt200_400;
   }
   if(dataset.Contains("Wqq_pt-400-600")){
      xs = XS_Wto2Q_pt400_600;
      xs_err = XSerr_Wto2Q_pt400_600;
   }
   if(dataset.Contains("Wqq_pt-600")){
      xs = XS_Wto2Q_pt600;
      xs_err = XSerr_Wto2Q_pt600;
   }
   if(dataset.Contains("VBFZqq")){
      xs = XS_VBFZto2Q;
      xs_err = XSerr_VBFZto2Q;
   }
   if(dataset.Contains("VBFWqq")){
      xs = XS_VBFWto2Q;
      xs_err = XSerr_VBFWto2Q;
   }
   return std::make_pair(xs, xs_err);

}

