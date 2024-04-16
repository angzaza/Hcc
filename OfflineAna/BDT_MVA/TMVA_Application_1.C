#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <algorithm>
#include "TString.h"

using namespace std;

void TMVA_Application_1(){

  TString fout_path = "TMVA_BDTG_2023C_out.root";
	TString method="BDTG";

  TString TMVA_weightfilename="dataset_v3Mar24/weights/TMVA_new_BDTG.weights.xml";

  float cutA=0.8;
  float cutB=0.9;
  float cutC=0.95;
  float cutD=1.0;

  double CvsB_jetC1, CvsB_jetC2, CvsL_jetC1, CvsL_jetC2, mqq, Deta_qq, Dphi_qq, Alpha_qq, pz_4jets, pt_norm, DR_HiggsVBF1, DR_HiggsVBF2, Dphi_qq_cc, jetEne_sum, jetPt_sum, QvsG_VBF1, QvsG_VBF2;
  int njets;
  double mCC;
  double mCC_BDTLoose, mCC_BDTMedium, mCC_BDTTight, mass_cJets, BDT_score;

  TMVA::Tools::Instance();
  TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
  float  forBDTevaluation1; 
  float  forBDTevaluation2; 
  float  forBDTevaluation3; 
  float  forBDTevaluation4; 
  float  forBDTevaluation5; 
  float  forBDTevaluation6; 
  float  forBDTevaluation7;  
  float  forBDTevaluation8;  
  float  forBDTevaluation9;  
  float  forBDTevaluation10;  
  float  forBDTevaluation11;  
  float  forBDTevaluation12;  
  float  forBDTevaluation13;  
  float  forBDTevaluation14;  
  float  forBDTevaluation15;  
  float  forBDTevaluation16;
  float  forBDTevaluation17;
  float  forBDTevaluation18;

  reader->TMVA::Reader::AddVariable( "CvsB_jetC1", &forBDTevaluation1 );
  reader->TMVA::Reader::AddVariable( "CvsB_jetC2", &forBDTevaluation2 );
  reader->TMVA::Reader::AddVariable( "CvsL_jetC1", &forBDTevaluation3 );
  reader->TMVA::Reader::AddVariable( "CvsL_jetC2", &forBDTevaluation4 );
  reader->TMVA::Reader::AddVariable( "mqq", &forBDTevaluation5 );
  reader->TMVA::Reader::AddVariable( "Deta_qq", &forBDTevaluation6 );
  reader->TMVA::Reader::AddVariable( "Dphi_qq", &forBDTevaluation7 );
  reader->TMVA::Reader::AddVariable( "Alpha_qq", &forBDTevaluation8 );
  reader->TMVA::Reader::AddVariable("QvsG_VBF1", &forBDTevaluation9);
  reader->TMVA::Reader::AddVariable("QvsG_VBF2", &forBDTevaluation10);
  reader->TMVA::Reader::AddVariable( "pz_4jets", &forBDTevaluation11 );
  reader->TMVA::Reader::AddVariable( "pt_norm", &forBDTevaluation12 );
  reader->TMVA::Reader::AddVariable( "DR_HiggsVBF1", &forBDTevaluation13 );
  reader->TMVA::Reader::AddVariable( "DR_HiggsVBF2", &forBDTevaluation14 );
  reader->TMVA::Reader::AddVariable( "Dphi_qq_cc", &forBDTevaluation15 );
  reader->TMVA::Reader::AddVariable( "njets", &forBDTevaluation16 );
  reader->TMVA::Reader::AddVariable( "jetEne_sum", &forBDTevaluation17 );
  reader->TMVA::Reader::AddVariable( "jetPt_sum", &forBDTevaluation18 );
  reader->TMVA::Reader::BookMVA("BDTG", TMVA_weightfilename);
  cout<<"Using weights in "<<TMVA_weightfilename<<endl;


  //TFile *f_bkg = new TFile("/lustrehome/azaza/HccAnalysis/CMSSW_12_4_3/src/Analysis/2023C_Hcc_ReReco22Sep_v3Mar2024/AnalysedTree_forAnalysis_data_2023C_Hcc.root");
  //TFile *f_bkg = new TFile("/lustrehome/azaza/HccAnalysis/CMSSW_12_4_3/src/Analysis/mergedFiles/Hcc_v3Mar24/AnalysedTree_data_2023C_Hcc_merged.root");
  TFile *f_bkg = new TFile("/lustrehome/azaza/HccAnalysis/CMSSW_12_4_3/src/Analysis/mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_VBFHCC_Hcc_merged.root");
  //TFile *f_bkg = new TFile("AnalysedTree_forAnalysis_data_2023C_Hcc_v3.root");
  TTree *t_bkg = (TTree*)f_bkg->Get("FinalTree");
  
  t_bkg->SetBranchAddress("CvsB_jetC1",&CvsB_jetC1);
  t_bkg->SetBranchAddress("CvsB_jetC2",&CvsB_jetC2);
  t_bkg->SetBranchAddress("CvsL_jetC1",&CvsL_jetC1);
  t_bkg->SetBranchAddress("CvsL_jetC2",&CvsL_jetC2);
  t_bkg->SetBranchAddress("mqq",&mqq);
  t_bkg->SetBranchAddress("Deta_qq",&Deta_qq);
  t_bkg->SetBranchAddress("Dphi_qq",&Dphi_qq);
  t_bkg->SetBranchAddress("Alpha_qq",&Alpha_qq);
  t_bkg->SetBranchAddress("QvsG_VBF1",&QvsG_VBF1);
  t_bkg->SetBranchAddress("QvsG_VBF2",&QvsG_VBF2);
  t_bkg->SetBranchAddress("pz_4jets",&pz_4jets);
  t_bkg->SetBranchAddress("pt_norm", &pt_norm);
  t_bkg->SetBranchAddress("DR_HiggsVBF1", &DR_HiggsVBF1);
  t_bkg->SetBranchAddress("DR_HiggsVBF2", &DR_HiggsVBF2);
  t_bkg->SetBranchAddress("Dphi_qq_cc", &Dphi_qq_cc);
  t_bkg->SetBranchAddress("njets", &njets);
  t_bkg->SetBranchAddress("jetEne_sum", &jetEne_sum);
  t_bkg->SetBranchAddress("jetPt_sum", &jetPt_sum);
  t_bkg->SetBranchAddress("mCC", &mCC);
  //t_bkg->SetBranchAddress("XS", &XS);
  //t_bkg->SetBranchAddress("Lumi", &Lumi);

  TString filename = "InvariantMass_signal_prova.root";
  TFile *hfile =0;
  hfile= TFile::Open(filename,"RECREATE");
  TTree *tree = new TTree("T","InvariantMass");
  tree->Branch("mass_cJets",&mass_cJets,"mass_cJets/D"); 
  tree->Branch("BDT_score",&BDT_score,"BDT_score/D");
  //tree->Branch("mCC_BDTLoose",&mCC_BDTLoose,"mCC_BDTLoose/D"); 
  //tree->Branch("mCC_BDTMedium",&mCC_BDTMedium,"mCC_BDTMedium/D"); 
  //tree->Branch("mCC_BDTTight",&mCC_BDTMedium,"mCC_BDTMedium/D"); 
  
  for(Long64_t ievt=0; ievt<t_bkg->GetEntries();ievt++){ 
 			t_bkg->GetEntry(ievt);
      forBDTevaluation1= CvsB_jetC1;
        forBDTevaluation2= CvsB_jetC2;
        forBDTevaluation3= CvsL_jetC1;
        forBDTevaluation4= CvsL_jetC2;
        forBDTevaluation5= mqq;
        forBDTevaluation6= Deta_qq;
        forBDTevaluation7= Dphi_qq;
        forBDTevaluation8= Alpha_qq;
        forBDTevaluation9= QvsG_VBF1;
        forBDTevaluation10= QvsG_VBF2;
        forBDTevaluation11= pz_4jets;
        forBDTevaluation12= pt_norm;
        forBDTevaluation13= DR_HiggsVBF1;
        forBDTevaluation14= DR_HiggsVBF2;
        forBDTevaluation15= Dphi_qq_cc;
        forBDTevaluation16= njets;
        forBDTevaluation17= jetEne_sum;
        forBDTevaluation18= jetPt_sum;

      mCC_BDTLoose=-1.; mCC_BDTMedium=-1.; mCC_BDTTight=-1.;

      double BDTG_decision = reader->EvaluateMVA(method);
      mass_cJets=mCC;
      BDT_score=BDTG_decision;
      //cout<<"BDT decision: "<<BDTG_decision<<endl;
      /*if(BDTG_decision>=cutA && BDTG_decision<cutB){ 
				mCC_BDTLoose=mCC;
        //cout<<"mass: "<<mCC<<endl;
      }
      if(BDTG_decision>=cutB && BDTG_decision<cutC) mCC_BDTMedium=mCC;
      if(BDTG_decision>=cutC && BDTG_decision<cutD) mCC_BDTTight=mCC;*/
      tree->Fill();
  }
      

  tree->Print();
  tree->Write();


}
