#define myAnalyzer_cxx
#define NCUTS 6

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TMVA/Reader.h>


void myAnalizer::TreeFin_Init(TTree *&tree, Double_t &isMC, Double_t &lumi, Double_t &run, Double_t &evt,int &entry, Double_t &puFactor, Double_t &pt_jetC1, Double_t &pt_jetC2, Double_t &pt_jetVBF1, Double_t &pt_jetVBF2, Double_t &eta_jetC1, Double_t &eta_jetC2, Double_t &eta_jetVBF1, Double_t &eta_jetVBF2, Double_t &CvsAll_jetC1, Double_t &CvsAll_jetC2, Double_t &CvsB_jetC1, Double_t &CvsB_jetC2, Double_t &CvsL_jetC1, Double_t &CvsL_jetC2, Double_t &mqq, Double_t &Deta_qq, Double_t &Dphi_qq, Double_t  &Alpha_qq, Double_t &qgl_VBF1, Double_t &qgl_VBF2, Double_t &pz_4jets, Double_t &pt_norm, Double_t &DR_HiggsVBF1, Double_t &DR_HiggsVBF2, Double_t &Dphi_qq_cc, int &njets, Double_t &jetEne_sum, Double_t  &jetPt_sum, Double_t &mCC){
		// Set tree branches
		tree->Branch("isMC", &isMC);
		tree->Branch("lumi", &lumi);
		tree->Branch("run", &run);
		tree->Branch("evt", &evt);
		tree->Branch("entry", &entry);
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
    listCut[2] = "pT_cut";
    listCut[3] = "2C_cands";
    listCut[4] = "CvsAll_0p5";
    listCut[5] = "VBF_cut";
  
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

