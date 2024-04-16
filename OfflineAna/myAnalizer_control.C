#define myAnalizer_control_cxx
#include "myAnalizer_control.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>

using namespace std;

struct jet_struct_2{
  TLorentzVector jet_p4;
  float CvsAll;
 // float CvsL;
 // float CvsB;
  bool c_tagged;
};

bool compareByCvsAll_2(const jet_struct_2 &a, const jet_struct_2 &b){

  if (std::isnan(a.CvsAll)) {
    // If a is NaN, b should come before a
    return false;
  }
  else if (std::isnan(b.CvsAll)) {
    // If b is NaN, a should come before b
    return true;
  }
  else{
    // Neither a nor b is NaN, compare their ctag values
    return a.CvsAll > b.CvsAll;
  }
}

bool ApplyVeto_SF(TH2D* histo, double eta, double phi){
  //cout<<"apply veto function aa"<<endl; 
  int binX = histo->GetXaxis()->FindBin(eta);
  int binY = histo->GetYaxis()->FindBin(phi);
  //cout<<"apply veto function bb"<<endl;
  if (histo->GetBinContent(binX, binY) != 0) return false;
  //cout<<"apply veto function cc"<<endl;
  return true;
}
float Get_PUweight_SF(TH1F* histo, int pu){
  int binX = histo->GetXaxis()->FindBin(pu);
  float puw = histo->GetBinContent(binX);
  //cout<<"apply veto function cc"<<endl;
  return puw;
}

void myAnalizer_control::Loop_control(TString type, TString datasetName, TString era)
{
//   In a ROOT session, you can do:
//      root> .L myAnalizer_control.C
//      root> myAnalizer_control t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
  
   //output file definition
   TString root_fileName = fileName;
   TFile *fout = new TFile(root_fileName, "RECREATE");
   fout->cd();
   TTree *tree = new TTree("FinalTree","FinalTree");
   cout<<"fileName: "<<root_fileName<<endl; 
   
  //name of the jet veto map file
  TString vetoMapFile;
  TFile* vetoFile;
  TH2D* vetoMapHisto;

  TFile* PUrewFile;
  TH1F* PUrewHisto;

   //variables
   bool verbose = false;
   bool Apply_jetVeto=true;
   bool Apply_PUrew=true;
   bool HLT_passed=false;
   bool HLT_PFJet80_passed=false;
   bool HLT_PFJet60_passed=false;
   bool HLT_Signal=false;
   bool HLT_Control=false;
   bool tag_passed=false;
   bool tag_matched=false;
   bool probe1_passed=false, probe2_passed=false, probe3_passed=false;
   TVector3 AK4puppi_jet1, AK4puppi_jet2, AK4puppi_jet3;
   TVector3 hltjet_p3;
   std::vector<jet_struct_2> AK4puppi_PNet;
   std::vector<jet_struct_2> AK4puppi_sel;
   std::vector<jet_struct_2> AK4puppi_cHiggs;
   std::vector<jet_struct_2> AK4puppi_VBF;
   TVector3 HLTjet80_p3, HLTjet60_p3, HLTAK4PFJetsLooseID_p3;
   double DeltaR_tag;
   double DeltaR_probe;

   double pt_bins[43];
   for(unsigned int i=0;i<37;i++){
     pt_bins[i]=i*5;
   }
   pt_bins[37]=200; pt_bins[38]=240; pt_bins[39]=280; pt_bins[40]=320; pt_bins[41]=360; pt_bins[42]=400;

   /*double pt_bins[20];
    for(unsigned int i=0;i<14;i++){
      pt_bins[i]=40+i*10;
   }
   pt_bins[14]=200; pt_bins[15]=240; pt_bins[16]=280; pt_bins[17]=320; pt_bins[18]=360; pt_bins[19]=400;*/

   double eta_bins[5];
   eta_bins[0]=0.; eta_bins[1]=1.4; eta_bins[2]=2.4; eta_bins[3]=3.0; eta_bins[4]=4.7;
   //histos
   /*TH1F *histoPt_jet1_tag = new TH1F("pt_jet1_tag","pt_jet1_tag",42, pt_bins);
   TH1F *histoPt_jet2_tag = new TH1F("pt_jet2_tag","pt_jet2_tag",42, pt_bins);
   TH1F *histoPt_jet3_tag = new TH1F("pt_jet3_tag","pt_jet3_tag",42, pt_bins);
   TH1F *histoPt_jet1_probe = new TH1F("pt_jet1_probe","pt_jet1_probe",42, pt_bins);
   TH1F *histoPt_jet2_probe = new TH1F("pt_jet2_probe","pt_jet2_probe",42, pt_bins);
   TH1F *histoPt_jet3_probe = new TH1F("pt_jet3_probe","pt_jet3_probe",42, pt_bins);*/

   TH1F *histoPNet_signal = new TH1F("PNet_signal", "PNet_signal", 20, 0, 1);
   TH1F *histoPNet_control = new TH1F("PNet_control", "PNet_control", 20, 0, 1);
   histoPNet_signal->Sumw2();
   histoPNet_control->Sumw2();


   TH2F *histoPtEta_jet1_tag = new TH2F("pt_jet1_tag","pt_jet1_tag",42, pt_bins,4, eta_bins);
   histoPtEta_jet1_tag->Sumw2();
   TH2F *histoPtEta_jet2_tag = new TH2F("pt_jet2_tag","pt_jet2_tag",42, pt_bins,4, eta_bins);
   histoPtEta_jet2_tag->Sumw2();
   TH2F *histoPtEta_jet3_tag = new TH2F("pt_jet3_tag","pt_jet3_tag",42, pt_bins,4, eta_bins);
   histoPtEta_jet3_tag->Sumw2();
   TH2F *histoPtEta_jet1_probe = new TH2F("pt_jet1_probe","pt_jet1_probe",42, pt_bins,4, eta_bins);
   histoPtEta_jet1_probe->Sumw2();
   TH2F *histoPtEta_jet2_probe = new TH2F("pt_jet2_probe","pt_jet2_probe",42, pt_bins,4, eta_bins);
   histoPtEta_jet2_probe->Sumw2();
   TH2F *histoPtEta_jet3_probe = new TH2F("pt_jet3_probe","pt_jet3_probe",42, pt_bins,4, eta_bins);
   histoPtEta_jet3_probe->Sumw2();

   double isMC = -99;
   if(datasetName.Contains("2023")) isMC=0;


   if(era.Contains("2023C")) vetoMapFile="/lustrehome/azaza/HccAnalysis/CMSSW_12_4_3/src/Analysis/jet_VetoMaps/Summer23Prompt23_RunC_v1.root";
   if(Apply_jetVeto==true){
     vetoFile = new TFile(vetoMapFile);
     if (!vetoFile || vetoFile->IsZombie()) {
       std::cerr << "Error: Could not open veto map file!" << std::endl;
       return;
     }

     vetoMapHisto = dynamic_cast<TH2D*>(vetoFile->Get("jetvetomap"));
     if (!vetoMapHisto) {
       std::cerr << "Error: Could not retrieve TH2D histogram from file!" << std::endl;
       vetoFile->Close();
       return;
      }
   }

   if(isMC!=0 && Apply_PUrew==true){
     PUrewFile = new TFile("/lustrehome/azaza/HccAnalysis/CMSSW_12_4_3/src/Analysis/PU_weights.root");
     if (!PUrewFile || PUrewFile->IsZombie()) {
       std::cerr << "Error: Could not open veto map file!" << std::endl;
       return;
     }

     PUrewHisto = dynamic_cast<TH1F*>(PUrewFile->Get("weights"));
     if (!PUrewHisto) {
       std::cerr << "Error: Could not retrieve TH1F histogram from file!" << std::endl;
       PUrewFile->Close();
       return;
      }
   }

   if (fChain == 0) return;

   //Long64_t nentries = 3000; // TO BE CHANGED
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if(verbose) cout<<"event:"<<jentry<<endl;
      if(verbose) cout<<"size hlt jets: "<<HLTAK4PFJetLoose_pt->size()<<endl;
      if(verbose) cout<<"size hltphjet60 jets: "<<HLTJet60_MatchedCalo_pt->size()<<endl;
      if(verbose) cout<<"size hltphjet80 jets: "<<HLTJet80_MatchedCalo_pt->size()<<endl;
      if(verbose) cout<<"size puppi jets: "<<AK4PuppiJets_pt->size()<<endl;
      if(verbose) cout<<"size PNet score: "<<jet_pfParticleNetAK4JetTags_probc->size()<<endl;
      HLT_passed=false;
      HLT_PFJet80_passed=false;
      HLT_PFJet60_passed=false;
      HLT_Signal=false;
      HLT_Control=false; 
      tag_passed=false;
      tag_matched=false;
      probe1_passed=false; probe2_passed=false; probe3_passed=false;
      AK4puppi_PNet.clear();      
      AK4puppi_sel.clear();      
      AK4puppi_cHiggs.clear();      
      AK4puppi_VBF.clear();      
      //HLT trigger
      for(int h=0; h<Trigger_hltname->size(); h++){
        TString hltName = Trigger_hltname->at(h);
        //if(Trigger_hltdecision->at(h)==1){
          //cout<<"hlt path: "<< hltName<<endl;
        //}
        if(hltName.Contains("HLT_PFJet80_v")   && Trigger_hltdecision->at(h) == 1){
          HLT_PFJet80_passed = true;
          HLT_passed=true;
          if(verbose) cout<<"Passed HLT_PFJet80" <<endl;
        }
        if(hltName.Contains("HLT_PFJet60_v")  && Trigger_hltdecision->at(h) == 1){
          HLT_PFJet60_passed = true;
          HLT_passed=true;
          if(verbose) cout<<"Passed HLT_PFJet60" <<endl;
        }
        if(hltName.Contains("HLT_QuadPFJet100_88_70_30_PNet1CvsAll0p5_VBF3Tight_v")  && Trigger_hltdecision->at(h) == 1){
          HLT_Signal = true;
          if(verbose) cout<<"Passed HLT_signal" <<endl;
        }
        if(hltName.Contains("HLT_QuadPFJet100_88_70_30_v")  && Trigger_hltdecision->at(h) == 1){
          HLT_Control = true;
          if(verbose) cout<<"Passed HLT_control" <<endl;
        }
      }

      //Pileup reweighting      
      float pu_weight=1.0;
      if(isMC!=0 && Apply_PUrew==true) pu_weight=Get_PUweight_SF(PUrewHisto, nInt);
     
      if(verbose) cout<<"applying jet veto maps" <<endl;
      //jet veto maps
      bool JetVetoed=false;
      for(unsigned int ijet=0; ijet< AK4PuppiJets_pt->size();ijet++){
            TVector3 AK4puppi_p3;
            AK4puppi_p3.SetPtEtaPhi(AK4PuppiJets_pt->at(ijet), AK4PuppiJets_eta->at(ijet), AK4PuppiJets_phi->at(ijet));
            //check if the jet contains a muon in dR<0.2
            float dRjmu_max=0.2;
            float dRjmu=-1;
            bool match_jmu=false;
            if(verbose) cout<<"matching with muons"<<endl;
            for(unsigned int j=0;j<Muon_pt->size();j++){
               TVector3 muon_p3;
               muon_p3.SetPtEtaPhi(Muon_pt->at(j),Muon_eta->at(j), Muon_phi->at(j));
               dRjmu = AK4puppi_p3.DeltaR(muon_p3);
               if(verbose) cout<<"dR jmu: "<<dRjmu<<endl;
               if(dRjmu<=dRjmu_max){
                 match_jmu=true;
                 if(verbose) cout<<"matched -- skipping current jet"<<endl;
                 break;
               }
            }
            if(match_jmu==true) continue;
            //cout<<"aa"<<endl;
            if(Apply_jetVeto==true){
              if(!ApplyVeto_SF(vetoMapHisto,AK4PuppiJets_eta->at(ijet), AK4PuppiJets_phi->at(ijet) )){
                //cout<<"aaaaa"<<endl; 
                JetVetoed=true;
                if(verbose) cout<<"event vetoed by jet maps"<<endl;
                break;
              }
            }
      }
      if(verbose) cout<<"finished applying jet veto maps" <<endl;
      if(JetVetoed==true) continue;

      if(verbose) cout<<"#### PNET SFs #############"<<endl;
      // PNet offline scale factors
      if(HLT_Control==true && AK4PuppiJets_pt->size()>=4){
         //bool match_hltjet=false;
         float pt_offlineThr[4]={105.,90.,75.,35.};
         float pt_hltThr[3]={100.,88.,70.};
         float dR_match[4]={-999.};
         bool match_hltjet[4]={false,false,false,true}; // the last is set to true because the match with hlt is not required for the last offline jet
      //if(HLT_Control==true && jet_pfParticleNetAK4JetTags_probc->size()>=4){  //to be subtituted by upper line once new ntuples are produced
        //cout<<"fired control path:"<<endl;
      //if((HLT_Signal == true || HLT_Control==true) && AK4PuppiJets_pt->size()>=4){

        for(unsigned int i=0; i<4; i++){
          if(verbose) cout<<"jet "<<i<<endl; 
          TLorentzVector AK4puppi_p4;
          TVector3 AK4puppi_p3;
          AK4puppi_p3.SetPtEtaPhi(AK4PuppiJets_pt->at(i), AK4PuppiJets_eta->at(i), AK4PuppiJets_phi->at(i));
          //check if the jet contains a muon in dR<0.2
          float dRjmu_max=0.2;
          float dRjmu=-1;
          bool match_jmu=false;
          if(verbose) cout<<"matching with muons"<<endl;
          for(unsigned int j=0;j<Muon_pt->size();j++){
            TVector3 muon_p3;
            muon_p3.SetPtEtaPhi(Muon_pt->at(j),Muon_eta->at(j), Muon_phi->at(j));
            dRjmu = AK4puppi_p3.DeltaR(muon_p3);
            if(verbose) cout<<"dR jmu: "<<dRjmu<<endl;
            if(dRjmu<=dRjmu_max){
              match_jmu=true;
              if(verbose) cout<<"matched -- skipping current jet"<<endl;
              break;
            }
          }
          if(match_jmu==true) break;
          if(verbose) cout<<"pt: "<<AK4PuppiJets_pt->at(i)<<"    eta: "<<AK4PuppiJets_eta->at(i)<<"      PNet: "<<jet_pfParticleNetAK4JetTags_probc->at(i)/(jet_pfParticleNetAK4JetTags_probc->at(i)+ jet_pfParticleNetAK4JetTags_probb->at(i) + jet_pfParticleNetAK4JetTags_probuds->at(i) + jet_pfParticleNetAK4JetTags_probg->at(i))<<endl;
          if(AK4PuppiJets_pt->at(i)<pt_offlineThr[i]) break; //check pt threshold offline jets
          int match_idx=-1;
          if(i==3) match_idx=1;
          float dR_max=0.4;
          if(i!=3){
            for(unsigned int j=0; j<HLTAK4PFJetLoose_pt->size(); j++){
              if(i==0 && verbose) cout<<"hlt jet "<<j<<"     pt: "<<HLTAK4PFJetLoose_pt->at(j)<<"      eta: "<<HLTAK4PFJetLoose_eta->at(j)<<"      phi: "<<HLTAK4PFJetLoose_phi->at(j)<<endl ;
              hltjet_p3.SetPtEtaPhi(HLTAK4PFJetLoose_pt->at(j), HLTAK4PFJetLoose_eta->at(j), HLTAK4PFJetLoose_phi->at(j));
              float this_dR=AK4puppi_p3.DeltaR(hltjet_p3);
              if(this_dR<dR_max){
                dR_max=this_dR;
                dR_match[i]=this_dR;
                //match_hltjet[i]=true;
                match_idx=j;
              }
            }

            if(match_idx>-1){
              if(verbose) cout<<"matched with hlt jet: "<<match_idx<<endl;
              hltjet_p3.SetPtEtaPhi(HLTAK4PFJetLoose_pt->at(match_idx), HLTAK4PFJetLoose_eta->at(match_idx), HLTAK4PFJetLoose_phi->at(match_idx));
              if(hltjet_p3.Pt()>pt_hltThr[i]) match_hltjet[i]=true;
              else break;
              if(verbose && match_hltjet[i]==true) cout<<"hlt jet above thr: "<<endl;
            }
          } 

          //cout<<"pt: "<<AK4PuppiJets_pt->at(i)<<"   eta: "<<AK4PuppiJets_eta->at(i)<<"   phi: "<<AK4PuppiJets_phi->at(i)<<"    mass: "<<AK4PuppiJets_mass->at(i)<<endl;
          //cout<<"probc: "<< jet_pfParticleNetAK4JetTags_probc->at(i)<<"   probb: "<<jet_pfParticleNetAK4JetTags_probb->at(i)<< "   probuds: "<<jet_pfParticleNetAK4JetTags_probuds->at(i)<<"     probg: "<<jet_pfParticleNetAK4JetTags_probg->at(i)<<endl;
          if(match_hltjet[i]==true) {
            //cout<<"matched with hlt jet "<<match_idx<<"    dR: "<<dR_match[i]<<endl;
            AK4puppi_p4.SetPtEtaPhiM(AK4PuppiJets_pt->at(i), AK4PuppiJets_eta->at(i), AK4PuppiJets_phi->at(i), AK4PuppiJets_mass->at(i));
            float num= jet_pfParticleNetAK4JetTags_probc->at(i);
            float den= jet_pfParticleNetAK4JetTags_probc->at(i)+ jet_pfParticleNetAK4JetTags_probb->at(i) + jet_pfParticleNetAK4JetTags_probuds->at(i) + jet_pfParticleNetAK4JetTags_probg->at(i);
            float CvsAll_val= num/den;
            if(verbose) cout<<"CvsAll jet"<<i<<": "<<CvsAll_val<<endl;
            //float CvsAll_val= jet_pfParticleNetAK4JetTags_CvsAll;
            jet_struct_2 jet_i={AK4puppi_p4,CvsAll_val, false};
            if(fabs(AK4puppi_p4.Eta())<4.7) AK4puppi_PNet.push_back(jet_i);
          }
          //cout<<"bbb"<<endl;
        } //end cycle on offline jets
        
        if(AK4puppi_PNet.size()>=4){
          std::sort(AK4puppi_PNet.begin(),AK4puppi_PNet.end(),compareByCvsAll_2 );
          int n_cjets=0;
          for(unsigned int k=0; k<4; k++){
            if(n_cjets<2){ //DA CONTROLLARE
              if(fabs((AK4puppi_PNet.at(k).jet_p4).Eta())<=2.4){
                AK4puppi_PNet.at(k).c_tagged=true;
                n_cjets++;
              }
            }
          }
          if(n_cjets==2){
            if(verbose) cout<<"2 cjets"<<endl;
            for(unsigned int j=0; j<4; j++){
              if(AK4puppi_PNet.at(j).c_tagged==true){
                AK4puppi_cHiggs.push_back(AK4puppi_PNet.at(j));
              }
              else{
                AK4puppi_VBF.push_back(AK4puppi_PNet.at(j));
              }
            } 
          }

          //cout<<"size Higgs jets: "<<AK4puppi_cHiggs.size()<<"         size VBF jets: "<<AK4puppi_VBF.size()<<endl;
          //cout<<"mass VBF jets: "<<(AK4puppi_VBF[0].jet_p4+AK4puppi_VBF[1].jet_p4).M()<<"    dEta: "<<fabs(AK4puppi_VBF[0].jet_p4.Eta()-AK4puppi_VBF[1].jet_p4.Eta())<<endl;

          //if(AK4puppi_PNet.size()>0 && match_hltjet[0]==true && match_hltjet[1]==true && match_hltjet[2]==true && match_hltjet[3]==true){
          if(AK4puppi_cHiggs.size()==2 && AK4puppi_VBF.size()==2 && (AK4puppi_VBF[0].jet_p4+AK4puppi_VBF[1].jet_p4).M()>=500. && fabs(AK4puppi_VBF[0].jet_p4.Eta() -AK4puppi_VBF[1].jet_p4.Eta())>=3.8){
            if(verbose) cout<<"passed VBF selection: "<<endl;
            histoPNet_control->Fill(AK4puppi_cHiggs.at(0).CvsAll,pu_weight*genWeight);
            //cout<<"histogram filled with CvsAll: "<<AK4puppi_cHiggs.at(0).CvsAll<<endl;
            if(HLT_Signal==true){
              //cout<<"fired hlt signal"<<endl;
              histoPNet_signal->Fill(AK4puppi_cHiggs.at(0).CvsAll,pu_weight*genWeight);
            }
          }
        }



				} 
      

      if(verbose) cout<<"#### Pt SFs #############"<<endl;
      // pt scale factors
      if(verbose) cout<<"HLT_passed: "<<HLT_passed<<endl;
      if(HLT_passed==true && AK4PuppiJets_pt->size()>=3){
      //if(HLT_passed==true && jet_pfParticleNetAK4JetTags_probc->size()>=3){ //TO BE CHANGED WITH UPPER LINE
        // tag selection
        // pT-leading offline jet should match an HLT jet with pT>130 GeV and eta<2.2
        int count_jetvetoed = 0;
        if(verbose) cout<<"starting cycle on the first 3 offline jets"<<endl;
        bool muveto_evt = false;
        for(unsigned int i=0; i<3; i++){
          if(verbose) cout<<"jet "<<i<<endl;
          TLorentzVector AK4puppi_p4;
          TVector3 AK4puppi_p3;
          AK4puppi_p3.SetPtEtaPhi(AK4PuppiJets_pt->at(i), AK4PuppiJets_eta->at(i), AK4PuppiJets_phi->at(i));
         //check if the jet contains a muon in dR<0.2
         float dRjmu_max=0.2;
         float dRjmu=-1;
         bool match_jmu=false;
         if(verbose) cout<<"matching with muons"<<endl;
         for(unsigned int j=0;j<Muon_pt->size();j++){
           TVector3 muon_p3;
           muon_p3.SetPtEtaPhi(Muon_pt->at(j),Muon_eta->at(j), Muon_phi->at(j));
           dRjmu = AK4puppi_p3.DeltaR(muon_p3);
           if(verbose) cout<<"dR jmu: "<<dRjmu<<endl;
           if(dRjmu<=dRjmu_max){
             match_jmu=true;
             if(verbose) cout<<"matched -- skipping current jet"<<endl;
             break;
           }
         }

        if(match_jmu==true && i<2){
            break; //se uno dei primi 2 jet contiene un muone devo scartare l'evento
            muveto_evt=true;
        }
        if(verbose) cout<<"pt: "<<AK4PuppiJets_pt->at(i)<<"    eta: "<<AK4PuppiJets_eta->at(i)<<"   phi: "<<AK4PuppiJets_phi->at(i)<<endl;
        }
        if(muveto_evt==true) continue;
        //cout<<"cc"<<endl;
        AK4puppi_jet1.SetPtEtaPhi(AK4PuppiJets_pt->at(0), AK4PuppiJets_eta->at(0), AK4PuppiJets_phi->at(0)); //offline leading jet    ---> tag
        AK4puppi_jet2.SetPtEtaPhi(AK4PuppiJets_pt->at(1), AK4PuppiJets_eta->at(1), AK4PuppiJets_phi->at(1)); //offline subleading jet ---> probe    
        AK4puppi_jet3.SetPtEtaPhi(AK4PuppiJets_pt->at(2), AK4PuppiJets_eta->at(2), AK4PuppiJets_phi->at(2));     

        //cout<<"pt jet1: "<< AK4puppi_jet1.Pt()<<endl;   
        //cout<<"pt jet2: "<< AK4puppi_jet2.Pt()<<endl;   
        //cout<<"pt jet3: "<< AK4puppi_jet3.Pt()<<endl;   
  
        float DeltaPhi_TagProbe=AK4puppi_jet1.DeltaPhi(AK4puppi_jet2);
        float ptMean_TagProbe=(AK4puppi_jet1.Pt()+AK4puppi_jet2.Pt())/2;        

        float tag_pTthr=110.;
        if(HLT_PFJet60_passed==false && HLT_PFJet80_passed==true){
          tag_pTthr=130.;
        } 
        if(verbose) cout<<"checking tag jet thresholds:"<<endl;
        if(AK4puppi_jet1.Pt()>tag_pTthr && fabs(AK4puppi_jet1.Eta())<2.2){
          if(verbose) cout<<"tag jet passed thr"<<endl;
          //for(unsigned int i=0; i< HLTJet80_pt->size(); i++){
          if(HLT_PFJet60_passed==false && HLT_PFJet80_passed==true){
            if(verbose) cout<<"hltpfjet60 not passed"<<endl;
            for(unsigned int i=0; i<HLTJet80_MatchedCalo_pt->size(); i++){
              //HLTjet80_p3.SetPtEtaPhi(HLTJet80_pt->at(i), HLTJet80_eta->at(i),HLTJet80_phi->at(i));
              HLTjet80_p3.SetPtEtaPhi(HLTJet80_MatchedCalo_pt->at(i), HLTJet80_MatchedCalo_eta->at(i),HLTJet80_MatchedCalo_phi->at(i));
              if(verbose) cout<<"hlt 80   pt: "<<HLTJet80_MatchedCalo_pt->at(i)<<"   eta: "<<HLTJet80_MatchedCalo_eta->at(i)<<"  phi: "<<HLTJet80_MatchedCalo_phi->at(i)<<endl; 
              DeltaR_tag=AK4puppi_jet1.DeltaR(HLTjet80_p3);
              if(DeltaR_tag<0.4){
                tag_matched=true;
              }    
            }
          }
          if(HLT_PFJet60_passed==true){
            for(unsigned int i=0; i<HLTJet60_MatchedCalo_pt->size(); i++){
              //HLTjet80_p3.SetPtEtaPhi(HLTJet80_pt->at(i), HLTJet80_eta->at(i),HLTJet80_phi->at(i));
              HLTjet60_p3.SetPtEtaPhi(HLTJet60_MatchedCalo_pt->at(i), HLTJet60_MatchedCalo_eta->at(i),HLTJet60_MatchedCalo_phi->at(i));
              if(verbose) cout<<"hlt 60   pt: "<<HLTJet60_MatchedCalo_pt->at(i)<<"   eta: "<<HLTJet60_MatchedCalo_eta->at(i)<<"  phi: "<<HLTJet60_MatchedCalo_phi->at(i)<<endl; 
              DeltaR_tag=AK4puppi_jet1.DeltaR(HLTjet60_p3);
              if(DeltaR_tag<0.4){
                tag_matched=true;
              }    
            }
          }
        }
        if(tag_matched && verbose) cout<<"tag matched with hlt"<<endl;
        if(verbose) cout<<"dPhi: "<<DeltaPhi_TagProbe<<"  0.3*ptMean_TagProbe: "<<0.3*ptMean_TagProbe<<endl;
        if(tag_matched==true && DeltaPhi_TagProbe>2.5 && AK4puppi_jet3.Pt()<0.3*ptMean_TagProbe){
          tag_passed=true;
          if(verbose) cout<<"tag passed"<<endl;
        }

        if(tag_passed==true){
          float DeltaR_min=0.4;
          int idx_match=-1;

          //cout<<"pt jet1: "<< AK4puppi_jet1.Pt()<<endl;   
          //cout<<"pt jet2: "<< AK4puppi_jet2.Pt()<<endl;   
          //cout<<"pt jet3: "<< AK4puppi_jet3.Pt()<<endl;   
          histoPtEta_jet1_tag->Fill(AK4puppi_jet2.Pt(), fabs(AK4puppi_jet2.Eta()),pu_weight*genWeight); //denominator
          histoPtEta_jet2_tag->Fill(AK4puppi_jet2.Pt(), fabs(AK4puppi_jet2.Eta()),pu_weight*genWeight); //denominator
          if(HLT_PFJet60_passed==true){
            histoPtEta_jet3_tag->Fill(AK4puppi_jet2.Pt(), fabs(AK4puppi_jet2.Eta()),pu_weight*genWeight); //denominator
          }
          //cout<<"size hlt jets: "<<HLTAK4PFJetLoose_pt->size()<<endl;
          cout<<"starting probe matching"<<endl;
          if(HLT_PFJet60_passed==true){
            for(unsigned int i=0; i<HLTJet60_MatchedCalo_pt->size(); i++){  //collezione da modificare
              HLTAK4PFJetsLooseID_p3.SetPtEtaPhi(HLTJet60_MatchedCalo_pt->at(i), HLTJet60_MatchedCalo_eta->at(i), HLTJet60_MatchedCalo_phi->at(i));
              if(verbose) cout<<"hlt jet"<<i<<"  pt: "<<HLTAK4PFJetsLooseID_p3.Pt()<<"  eta: "<<HLTAK4PFJetsLooseID_p3.Eta()<<"  phi: "<< HLTAK4PFJetsLooseID_p3.Phi() <<endl;

              DeltaR_probe=AK4puppi_jet2.DeltaR(HLTAK4PFJetsLooseID_p3);
              if(verbose) cout<<"dR: "<<DeltaR_probe<<endl;
              if(DeltaR_probe<DeltaR_min){
                DeltaR_min=DeltaR_probe;
                idx_match=i;
              }
            }
          }

          if(HLT_PFJet60_passed==false && HLT_PFJet80_passed==true){
            for(unsigned int i=0; i<HLTJet80_MatchedCalo_pt->size(); i++){  //collezione da modificare
              HLTAK4PFJetsLooseID_p3.SetPtEtaPhi(HLTJet80_MatchedCalo_pt->at(i), HLTJet80_MatchedCalo_eta->at(i), HLTJet80_MatchedCalo_phi->at(i));
              if(verbose) cout<<"hlt jet"<<i<<"  pt: "<<HLTAK4PFJetsLooseID_p3.Pt()<<"  eta: "<<HLTAK4PFJetsLooseID_p3.Eta()<<"  phi: "<< HLTAK4PFJetsLooseID_p3.Phi() <<endl;
              //cout<<"pt hlt probe: "<<HLTAK4PFJetsLooseID_p3.Pt()<<endl;

              DeltaR_probe=AK4puppi_jet2.DeltaR(HLTAK4PFJetsLooseID_p3);
              if(verbose) cout<<"dR: "<<DeltaR_probe<<endl;
              if(DeltaR_probe<DeltaR_min){
                DeltaR_min=DeltaR_probe;
                idx_match=i;
              }
            }
          }
          if(idx_match!=-1){
              if(verbose) cout<<"probe matched with hlt jet"<<idx_match<<endl;
              if(HLT_PFJet60_passed==true){
                HLTAK4PFJetsLooseID_p3.SetPtEtaPhi(HLTJet60_MatchedCalo_pt->at(idx_match), HLTJet60_MatchedCalo_eta->at(idx_match), HLTJet60_MatchedCalo_phi->at(idx_match));
                if(HLTAK4PFJetsLooseID_p3.Pt()>70.) probe3_passed=true;
              }
              if(HLT_PFJet60_passed==false && HLT_PFJet80_passed==true){
                HLTAK4PFJetsLooseID_p3.SetPtEtaPhi(HLTJet80_MatchedCalo_pt->at(idx_match), HLTJet80_MatchedCalo_eta->at(idx_match), HLTJet80_MatchedCalo_phi->at(idx_match));
              }
              //if(HLTAK4PFJetsLooseID_p3.Pt()>70.) probe3_passed=true;
              if(HLTAK4PFJetsLooseID_p3.Pt()>88.) probe2_passed=true;
              if(HLTAK4PFJetsLooseID_p3.Pt()>100.) probe1_passed=true;


          }
          if(probe3_passed==true){ 
            histoPtEta_jet3_probe->Fill(AK4puppi_jet2.Pt(), fabs(AK4puppi_jet2.Eta()),pu_weight*genWeight);
           if(verbose) cout<<"probe 3 passed"<<endl;
          }
          if(probe2_passed==true){ 
            histoPtEta_jet2_probe->Fill(AK4puppi_jet2.Pt(), fabs(AK4puppi_jet2.Eta()),pu_weight*genWeight);
            if(verbose) cout<<"probe 2 passed"<<endl;
          } 
          if(probe1_passed==true){ 
            histoPtEta_jet1_probe->Fill(AK4puppi_jet2.Pt(), fabs(AK4puppi_jet2.Eta()),pu_weight*genWeight);
            if(verbose) cout<<"probe 1 passed"<<endl; 
         }
        }
         
      }
  
   } //end cycle on events

   if(Apply_jetVeto==true) vetoFile->Close();
   if(isMC!=0 && Apply_PUrew==true) PUrewFile->Close();

   fout->cd();
   TH1D* histoPt_jet1_tag[4];
   TH1D* histoPt_jet2_tag[4];
   TH1D* histoPt_jet3_tag[4];

   TH1D* histoPt_jet1_probe[4];
   TH1D* histoPt_jet2_probe[4];
   TH1D* histoPt_jet3_probe[4];

   TEfficiency* ptEff_jet1[4];
   TEfficiency* ptEff_jet2[4];
   TEfficiency* ptEff_jet3[4];

   for(int i=0;i<4;i++){
     histoPt_jet1_tag[i]=histoPtEta_jet1_tag->ProjectionX(Form("jet1_tag_%d",i),i+1,i+1);
     histoPt_jet2_tag[i]=histoPtEta_jet2_tag->ProjectionX(Form("jet2_tag_%d",i),i+1,i+1);
     histoPt_jet3_tag[i]=histoPtEta_jet3_tag->ProjectionX(Form("jet3_tag_%d",i),i+1,i+1);
     histoPt_jet1_probe[i]=histoPtEta_jet1_probe->ProjectionX(Form("jet1_probe_%d",i),i+1,i+1);
     histoPt_jet2_probe[i]=histoPtEta_jet2_probe->ProjectionX(Form("jet2_probe_%d",i),i+1,i+1);
     histoPt_jet3_probe[i]=histoPtEta_jet3_probe->ProjectionX(Form("jet3_probe_%d",i),i+1,i+1);

     if(TEfficiency::CheckConsistency(*histoPt_jet1_probe[i],*histoPt_jet1_tag[i])){
       ptEff_jet1[i] = new TEfficiency(*histoPt_jet1_probe[i],*histoPt_jet1_tag[i]);
     }
    
     if(TEfficiency::CheckConsistency(*histoPt_jet2_probe[i],*histoPt_jet2_tag[i])){
       ptEff_jet2[i] = new TEfficiency(*histoPt_jet2_probe[i],*histoPt_jet2_tag[i]);
     }
     
     if(TEfficiency::CheckConsistency(*histoPt_jet3_probe[i],*histoPt_jet3_tag[i])){
       ptEff_jet3[i] = new TEfficiency(*histoPt_jet3_probe[i],*histoPt_jet3_tag[i]);
     }

     ptEff_jet1[i]->Write();
     ptEff_jet2[i]->Write();
     ptEff_jet3[i]->Write();
   }
   
   TEfficiency* PNet_Eff = 0;
   if(TEfficiency::CheckConsistency(*histoPNet_signal,*histoPNet_control)){
     PNet_Eff = new TEfficiency(*histoPNet_signal,*histoPNet_control);
   }
   
   PNet_Eff->Write();

   /*TEfficiency* ptEff_jet1 = 0;
   if(TEfficiency::CheckConsistency(*histoPt_jet1_probe,*histoPt_jet1_tag)){
     ptEff_jet1 = new TEfficiency(*histoPt_jet1_probe,*histoPt_jet1_tag);
   }
   ptEff_jet1->SetTitle("ptEff_jet1");
   
   TEfficiency* ptEff_jet2 = 0;
   if(TEfficiency::CheckConsistency(*histoPt_jet2_probe,*histoPt_jet2_tag)){
     ptEff_jet2 = new TEfficiency(*histoPt_jet2_probe,*histoPt_jet2_tag);
   }
   ptEff_jet2->SetTitle("ptEff_jet2");

   TEfficiency* ptEff_jet3 = 0;
   if(TEfficiency::CheckConsistency(*histoPt_jet3_probe,*histoPt_jet3_tag)){
     ptEff_jet3 = new TEfficiency(*histoPt_jet3_probe,*histoPt_jet3_tag);
   }
   ptEff_jet3->SetTitle("ptEff_jet3");*/


    
    //ptEff_jet1->Write();
    //ptEff_jet2->Write();
    //ptEff_jet3->Write();
    fout->Write();
    fout->Close();

}
