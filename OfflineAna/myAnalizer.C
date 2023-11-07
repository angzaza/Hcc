#define myAnalizer_cxx
#define NCUTS 6
#include "myAnalizer.h"
#include "Utilities.C"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;

bool compareByPt(const TLorentzVector &a, const TLorentzVector &b){
  return a.Pt() > b.Pt();
}

struct jet_struct{
  TLorentzVector jet_p4;
  float CvsAll;
  float CvsL;
  float CvsB;
  float QGL;
  bool c_tagged;
};

bool compareByCvsAll(const jet_struct &a, const jet_struct &b){
  return a.CvsAll > b.CvsAll;
}

bool compareByCvsL(const jet_struct &a, const jet_struct &b){
  return a.CvsL > b.CvsL;
}



void myAnalizer::Loop_Hcc(TString type, TString datasetName)
{
//   In a ROOT session, you can do:
//      root> .L myAnalizer.C
//      root> myAnalizer t
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


   //variable and histo declaration
   TLorentzVector hltjet_p4;
   TLorentzVector GENjet_p4;
   TLorentzVector AK4puppi_p4;
   vector<TLorentzVector> AK4puppi_p4_vec;
   vector<TLorentzVector> cjets_hlt, VBFjets_hlt, fakejets_hlt;
   TVector3 hltjet_p3, quark_p3, GENjet_p3;
   int match_index=-10;
   int match_hltGen_index=-10;
   vector<bool> match_free, match_gen_free;
   bool match_gen_quark=false;
   bool match_hlt_gen=false;
   bool trigger_VBFPNet=false;

   vector<float> hltjet_pt;

   int hlt_count=0;
   int hltdouble_count=0;
   int hltOR_count=0;

   bool PtThr_passed=false;

   //define a vector of jet_struct
   std::vector<jet_struct> AK4puppi_sel;
   std::vector<jet_struct> AK4puppi_cHiggs;
   std::vector<jet_struct> AK4puppi_VBF;
   bool passed_selection=false;
   //

   float ParticleNet_CvsAll_value, ParticleNet_CvsB_value, ParticleNet_CvsL_value;
   TH1I *hCutEffEvt = new TH1I("CutEff_NEvents", "CutEff_NEvents", NCUTS, 0.5, (NCUTS+0.5));

   TH1F* PassedEvents= new TH1F("passed events","passed events",2,-0.5,1.5);
   //TH1F* CvsAll_ptSel= new TH1F("PNet CvsAll (pt selection only)","PNet CvsAll (py selection only)",20,0,1);
   TH1F *histo_pt_jet1=new TH1F("pt_jet1","pt_jet1",100,0,400);
   TH1F *histo_pt_jet2=new TH1F("pt_jet2","pt_jet2",100,0,400);
   TH1F *histo_pt_jet3=new TH1F("pt_jet3","pt_jet3",100,0,400);
   TH1F *histo_pt_jet4=new TH1F("pt_jet4","pt_jet4",100,0,400);
   TH1F *histo_eta_jet1=new TH1F("eta_jet1"," eta_jet1",50,-5,5);
   TH1F *histo_eta_jet2=new TH1F("eta_jet2"," eta_jet2",50,-5,5);
   TH1F *histo_eta_jet3=new TH1F("eta_jet3"," eta_jet3",50,-5,5);
   TH1F *histo_eta_jet4=new TH1F("eta_jet4"," eta_jet4",50,-5,5);
   TH1F *histo_pt_jetC1=new TH1F("pt_jetc1","pt_jetc1",100,0,400);
   TH1F *histo_pt_jetC2=new TH1F("pt_jetc2","pt_jetc2",100,0,400);
   TH1F *histo_pt_jetVBF1=new TH1F("pt_jetVBF1","pt_jetVBF1",100,0,400);
   TH1F *histo_pt_jetVBF2=new TH1F("pt_jetVBF2","pt_jetVBF2",100,0,400);
   TH1F *histo_eta_jetC1=new TH1F("eta_jetc1"," eta_jetc1",50,-5,5);
   TH1F *histo_eta_jetC2=new TH1F("eta_jetc2"," eta_jetc2",50,-5,5);
   TH1F *histo_eta_jetVBF1=new TH1F("eta_jetVBF1"," eta_jetVBF2",50,-5,5);
   TH1F *histo_eta_jetVBF2=new TH1F("eta_jetVBF2"," eta_jetVBF2",50,-5,5);


   TH1F* histo_deltaEta_cjets= new TH1F("DeltaEta_cjets","DelatEta_cjets",50,0,12);
   TH1F* histo_deltaEta_VBFjets= new TH1F("DeltaEta_VBFjets","DeltaEta_VBFjets",50,0,12);

   TH1F* histo_invMass_cjets= new TH1F("InvMass_cjets","InvMass_cjets",100,0,400);
   TH1F* histo_invMass_VBFjets= new TH1F("InvMass_VBFjets","InvMass_VBFjets",100,0,1500);

   TH1F* histo_CvsAll_c1= new TH1F("PNet_CvsAll_jetc1","PNet_CvsAll_jetc1",20,0,1);
   TH1F* histo_CvsAll_c2= new TH1F("PNet_CvsAll_jetc2","PNet_CvsAll_jetc2",20,0,1);
   TH1F* histo_CvsL_c1= new TH1F("PNet_CvsL_jetc1","PNet_CvsL_jetc1",20,0,1);
   TH1F* histo_CvsL_c2= new TH1F("PNet_CvsL_jetc2","PNet_CvsL_jetc2",20,0,1);
   TH1F* histo_CvsB_c1= new TH1F("PNet_CvsB_jetc1","PNet_CvsB_jetc1",20,0,1);
   TH1F* histo_CvsB_c2= new TH1F("PNet_CvsB_jetc2","PNet_CvsB_jetc2",20,0,1);

   double isMC = -99;
   double run_n = 0, lumi_n = 0, evt_n = 0, pileupFactor=0;
   bool HLT_passed=false;

   double pt_jetC1=0, pt_jetC2=0, pt_jetVBF1=0, pt_jetVBF2=0, eta_jetC1=0, eta_jetC2=0, eta_jetVBF1=0, eta_jetVBF2=0, CvsAll_jetC1=0, CvsAll_jetC2=0, CvsB_jetC1=0, CvsB_jetC2=0, CvsL_jetC1=0, CvsL_jetC2=0;
   double mqq=0, Deta_qq=0, Dphi_qq=0, Alpha_qq=0, qgl_VBF1=0, qgl_VBF2=0, pz_4jets=0, pt_norm=0, DR_HiggsVBF1=0, DR_HiggsVBF2=0,  Dphi_qq_cc=0, jetEne_sum=0, jetPt_sum=0, mCC=0;
   int njets=0;
   int entry=0;
   int cutevt[NCUTS] = {0};

   TString listCut[NCUTS];
   Fill_CutName(listCut);

   //output file definition
   TString root_fileName = fileName;
   TFile *fout = new TFile(root_fileName, "RECREATE");
   fout->cd();
   TTree *tree = new TTree("FinalTree","FinalTree");
   //initialize output tree
   TreeFin_Init(tree, isMC, lumi_n, run_n, evt_n, entry, pileupFactor, pt_jetC1, pt_jetC2, pt_jetVBF1, pt_jetVBF2, eta_jetC1, eta_jetC2, eta_jetVBF1, eta_jetVBF2, CvsAll_jetC1, CvsAll_jetC2, CvsB_jetC1, CvsB_jetC2, CvsL_jetC1, CvsL_jetC2, mqq, Deta_qq, Dphi_qq, Alpha_qq, qgl_VBF1, qgl_VBF2, pz_4jets, pt_norm, DR_HiggsVBF1, DR_HiggsVBF2  ,Dphi_qq_cc, njets, jetEne_sum, jetPt_sum, mCC);

   if(datasetName.Contains("2023")) isMC=0;
    else{
        if(datasetName.Contains("QCD")) isMC=1;
       // if(datasetName.Contains("Bp")) isMC=2;
       // if(datasetName.Contains("B0")) isMC=3;
    }
   
   if (fChain == 0) return;

   //Long64_t nentries = 1000;
   Long64_t nentries = fChain->GetEntriesFast();
   //Long64_t nentries = 1000;

   Long64_t nbytes = 0, nb = 0;
   
   //cycle on events
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      entry=jentry;
     // cout<<"evt: "<<entry<<endl;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      bool Eff_counter[NCUTS] = {false};
      AK4puppi_sel.clear();
      AK4puppi_cHiggs.clear();
      AK4puppi_VBF.clear();
      HLT_passed=false;
    
      //set default values
      PtThr_passed=false;
      passed_selection=false;

      Eff_counter[0] = true;

      run_n = Run; lumi_n = LumiSect; evt_n = Event; pileupFactor=nInt;  
    
      njets=0; jetEne_sum=0; jetPt_sum=0;
 
      //HLT trigger
      for(int h=0; h<Trigger_hltname->size(); h++){
        TString hltName = Trigger_hltname->at(h);
        //if( (hltName.Contains("HLT_QuadPFJet103_88_75_15_v") || hltName.Contains("HLT_QuadPFJet105_88_76_15_v") || hltName.Contains("HLT_QuadPFJet111_90_80_15_v") || hltName.Contains("HLT PFJet80_v") || hltName.Contains("HLT_QuadPFJet103_88_75_15_PFBTagDeepJet_1p3_VBF2_v") || hltName.Contains("HLT_QuadPFJet105_88_76_15_PFBTagDeepJet_1p3_VBF2_v") || hltName.Contains("HLT_QuadPFJet111_90_80_15_PFBTagDeepJet_1p3_VBF2_v") || hltName.Contains("HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepJet_1p3_7p7_VBF1_v") || hltName.Contains("HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepJet_1p3_7p7_VBF1_v") || hltName.Contains("HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepJet_1p3_7p7_VBF1_v")) && Trigger_hltdecision->at(h) == 1){
        if( hltName.Contains("HLT_QuadPFJet100_88_70_30_PNet1CvsAll0p5_VBF3Tight_v") && Trigger_hltdecision->at(h) == 1 ){
           HLT_passed = true;
         }
      }

      if(!HLT_passed){ // se non passa HLT, aggiorna il conteggio degli eventi e passa all'evento successivo
        cutevt[0]++;
        continue;
      }

      if(HLT_passed==true){
        Eff_counter[1] = true;

        bool veto=false;
        //veto on electrons and muons
        for(unsigned int iele=0; iele< Ele_pt->size(); iele++){
          if(Ele_pt->at(iele)<7.0 || fabs(Ele_eta->at(iele)>2.5) || fabs(Ele_dxy->at(iele))>0.05 || fabs(Ele_dz->at(iele))>0.2) continue;
          if(Ele_IsoCal->at(iele)<=0.4){
            veto=true;
            break;
          }
        }
        if(veto==false){
          for(unsigned int imu=0; imu< Muon_pt->size(); imu++){
            if(Muon_pt->at(imu)<5.0 || fabs(Muon_eta->at(imu)>2.5) || fabs(Muon_dxy->at(imu))>0.25 || fabs(Muon_dz->at(imu))>1.0) continue;
            if(Muon_PF_Iso_R04->at(imu)<=0.4){
              veto=true;
              break;
            }
          }
        } 
        if(veto==false){
          //Eff_counter[2] = true;
          //cycle on AK4 puppi jets
          for(unsigned int ijet=0; ijet< AK4PuppiJets_pt->size();ijet++){
            bool match_hltjet=false;
            TVector3 AK4puppi_p3;
            AK4puppi_p3.SetPtEtaPhi(AK4PuppiJets_pt->at(ijet), AK4PuppiJets_eta->at(ijet), AK4PuppiJets_phi->at(ijet));
            //match with L3 trigger objects
            float dR_max=0.4;
            for(unsigned int j=0; j<HLTAK4PFJetLoose_pt->size(); j++){
              TVector3 hltjet_p3;
              hltjet_p3.SetPtEtaPhi(HLTAK4PFJetLoose_pt->at(j), HLTAK4PFJetLoose_eta->at(j), HLTAK4PFJetLoose_phi->at(j));
              float this_dR=AK4puppi_p3.DeltaR(hltjet_p3);
              if(this_dR<dR_max){
                match_hltjet=true;
                break;
              }
            }       

            if(match_hltjet==true){
              AK4puppi_p4.SetPtEtaPhiM(AK4PuppiJets_pt->at(ijet), AK4PuppiJets_eta->at(ijet), AK4PuppiJets_phi->at(ijet), AK4PuppiJets_mass->at(ijet));
             // cout<<"jet "<<ijet<<endl;
              float num= jet_pfParticleNetAK4JetTags_probc->at(ijet);
              float den= jet_pfParticleNetAK4JetTags_probc->at(ijet)+ jet_pfParticleNetAK4JetTags_probb->at(ijet) + jet_pfParticleNetAK4JetTags_probuds->at(ijet) + jet_pfParticleNetAK4JetTags_probg->at(ijet);
              //float CvsAll_val= (den!=0 || (num/den)<999)? num/den : -1;
              float CvsAll_val=num/den;
              den = jet_pfParticleNetAK4JetTags_probc->at(ijet) + jet_pfParticleNetAK4JetTags_probuds->at(ijet) + jet_pfParticleNetAK4JetTags_probg->at(ijet);
              float CvsL_val= den!=0 ? num/den : -1;
              //float CvsL_val= num/den;
              den = jet_pfParticleNetAK4JetTags_probc->at(ijet) + jet_pfParticleNetAK4JetTags_probb->at(ijet);
              float CvsB_val= den!=0 ? num/den : -1;
              //float CvsB_val= num/den;
              num = jet_pfParticleNetAK4JetTags_probuds->at(ijet);  
              den = jet_pfParticleNetAK4JetTags_probg->at(ijet)+ jet_pfParticleNetAK4JetTags_probuds->at(ijet);  
              float qgl= den!=0 ? num/den : -1;
              /*cout<<"probc: "<<jet_pfParticleNetAK4JetTags_probc->at(ijet)<<endl;
              cout<<"probb: "<<jet_pfParticleNetAK4JetTags_probb->at(ijet)<<endl;
              cout<<"probuds: "<<jet_pfParticleNetAK4JetTags_probuds->at(ijet)<<endl;
              cout<<"probg: "<<jet_pfParticleNetAK4JetTags_probg->at(ijet)<<endl;
              cout<<"CvsAll: "<<CvsAll_val<<endl;
              cout<<"CvsL: "<<CvsL_val<<endl;
              cout<<"CvsB: "<<CvsB_val<<endl;
              cout<<"QvsG: "<<qgl<<endl;*/
              //float CvsAll_val=jet_pfParticleNetAK4JetTags_CvsAll->at(ijet);
              //float CvsB_val=jet_pfParticleNetAK4JetTags_CvsB->at(ijet);
              //float CvsL_val=jet_pfParticleNetAK4JetTags_CvsL->at(ijet);
              //float qgl=jet_pfParticleNetAK4JetTags_QvsG->at(ijet);

              if(fabs(AK4puppi_p4.Eta())<4.7){
                jet_struct jet_i={AK4puppi_p4,CvsAll_val,CvsL_val,CvsB_val,qgl,false};
                AK4puppi_sel.push_back(jet_i);
              }

              if(AK4puppi_p4.Pt()>20. && fabs(AK4puppi_p4.Eta()<2.4)){
                njets++;
              }

              if(ijet>3 && AK4puppi_p4.Pt()>30. && fabs(AK4puppi_p4.Eta()<2.4)){
                jetEne_sum+=AK4puppi_p4.E();
                jetPt_sum+=AK4puppi_p4.Pt();
              }
            }//end condition on trigger matching
          }//end cycle on puppi
          
          if(AK4puppi_sel.size()>=4){
            //resize the vector of selected jets, keeping only the first 4 of them (pt sorted)
            AK4puppi_sel.resize(4);
            float pt0= ((AK4puppi_sel.at(0)).jet_p4).Pt();
            float pt1= ((AK4puppi_sel.at(1)).jet_p4).Pt();
            float pt2= ((AK4puppi_sel.at(2)).jet_p4).Pt();
            float pt3= ((AK4puppi_sel.at(3)).jet_p4).Pt();
            if(pt0>=105.0 && pt1>= 90. && pt2>=75. && pt3>=35.){
              PtThr_passed=true;
              histo_pt_jet1->Fill(pt0);
              histo_pt_jet2->Fill(pt1);
              histo_pt_jet3->Fill(pt2);
              histo_pt_jet4->Fill(pt3);
              histo_eta_jet1->Fill(((AK4puppi_sel.at(0)).jet_p4).Eta());
              histo_eta_jet2->Fill(((AK4puppi_sel.at(1)).jet_p4).Eta());
              histo_eta_jet3->Fill(((AK4puppi_sel.at(2)).jet_p4).Eta());
              histo_eta_jet4->Fill(((AK4puppi_sel.at(3)).jet_p4).Eta());
            }
            if(PtThr_passed==true){
             Eff_counter[2]=true;
             //sort the AK4puppi_sel by descending CvsAll score
             std::sort(AK4puppi_sel.begin(),AK4puppi_sel.end(),compareByCvsL );
             // CvsAll_ptSel->Fill(AK4puppi_sel.at(0).CvsAll);
             //tag as c jets the two jets with the highest CvsAll score within eta<2.4
             int n_cjets=0;
             for(unsigned int k=0; k<4; k++){
               if(n_cjets<2){ //DA CONTROLLARE
                 if(fabs((AK4puppi_sel.at(k).jet_p4).Eta())<=2.4){
                   AK4puppi_sel.at(k).c_tagged=true;
                   n_cjets++;
                 }
               }  
             }
             if(n_cjets==2){
               Eff_counter[3]=true;
               //cout<<"2 c-jets found"<<endl;
               for(unsigned int j=0; j<4; j++){
                 if(AK4puppi_sel.at(j).c_tagged==true){
                   AK4puppi_cHiggs.push_back(AK4puppi_sel.at(j));
                   /*cout<<"c jet:"<<endl;
                   cout<<"prob CvsAll:"<<AK4puppi_sel.at(j).CvsAll<<endl;
                   cout<<"prob CvsB:"<<AK4puppi_sel.at(j).CvsB<<endl;
                   cout<<"prob CvsL:"<<AK4puppi_sel.at(j).CvsL<<endl;
                   cout<<"prob QvsG:"<<AK4puppi_sel.at(j).QGL<<endl;*/
                 }
                 else{
                   AK4puppi_VBF.push_back(AK4puppi_sel.at(j));
                  /* cout<<"VBF jet"<<endl;
                   cout<<"prob CvsAll:"<<AK4puppi_sel.at(j).CvsAll<<endl;
                   cout<<"prob CvsB:"<<AK4puppi_sel.at(j).CvsB<<endl;
                   cout<<"prob CvsL:"<<AK4puppi_sel.at(j).CvsL<<endl;
                   cout<<"prob QvsG:"<<AK4puppi_sel.at(j).QGL<<endl;*/
                 }
               }
               //cout<<"CvsAll c1:"<<AK4puppi_cHiggs.at(0).CvsAll<<endl;
               if(AK4puppi_cHiggs.size()==2 && AK4puppi_VBF.size()==2){
                 //if(AK4puppi_cHiggs.at(0).CvsAll>=0.5){ //condition on ctag verified
                 Eff_counter[4]=true;
                 TLorentzVector jet_VBF1=AK4puppi_VBF.at(0).jet_p4;
                 TLorentzVector jet_VBF2=AK4puppi_VBF.at(1).jet_p4;
                 if((jet_VBF1+jet_VBF2).M()>=500. && fabs(jet_VBF1.Eta()-jet_VBF2.Eta())>=3.8){
                   passed_selection=true;
                 }
               //}
               }
             } 
           }
         }
       }

       //PassedEvents->Fill(passed_selection);

       if(passed_selection==true){
         Eff_counter[5]=true;
         TLorentzVector jet_VBF1=AK4puppi_VBF.at(0).jet_p4;
         TLorentzVector jet_VBF2=AK4puppi_VBF.at(1).jet_p4;
         TLorentzVector jet_c1=AK4puppi_cHiggs.at(0).jet_p4;
         TLorentzVector jet_c2=AK4puppi_cHiggs.at(1).jet_p4;
         pt_jetC1=jet_c1.Pt();
         pt_jetC2=jet_c2.Pt();
     
         pt_jetVBF1=jet_VBF1.Pt();
         pt_jetVBF2=jet_VBF2.Pt();

         eta_jetC1=jet_c1.Eta();
         eta_jetC2=jet_c2.Eta();

         eta_jetVBF1=jet_VBF1.Eta();
         eta_jetVBF2=jet_VBF2.Eta();

         histo_pt_jetC1->Fill(jet_c1.Pt());
         histo_pt_jetC2->Fill(jet_c2.Pt());

         histo_eta_jetC1->Fill(jet_c1.Eta());
         histo_eta_jetC2->Fill(jet_c2.Eta());

         histo_pt_jetVBF1->Fill(jet_VBF1.Pt());
         histo_pt_jetVBF2->Fill(jet_VBF2.Pt());

         histo_eta_jetVBF1->Fill(jet_VBF1.Eta());
         histo_eta_jetVBF2->Fill(jet_VBF2.Eta());
         histo_deltaEta_cjets->Fill(fabs(jet_c1.Eta()-jet_c2.Eta()));
         histo_deltaEta_VBFjets->Fill(fabs(jet_VBF1.Eta()-jet_VBF2.Eta()));

         histo_invMass_cjets->Fill((jet_c1+jet_c2).M());
         histo_invMass_VBFjets->Fill((jet_VBF1+jet_VBF2).M());

         histo_CvsAll_c1->Fill(AK4puppi_cHiggs.at(0).CvsAll);
         histo_CvsAll_c2->Fill(AK4puppi_cHiggs.at(1).CvsAll);

         histo_CvsL_c1->Fill(AK4puppi_cHiggs.at(0).CvsL);
         histo_CvsL_c2->Fill(AK4puppi_cHiggs.at(1).CvsL);

         histo_CvsB_c1->Fill(AK4puppi_cHiggs.at(0).CvsB);
         histo_CvsB_c2->Fill(AK4puppi_cHiggs.at(1).CvsB);

         // BDT variables
         //Deta_qq=0, Dphi_qq=0, Alfa_qq=0, qgl_qq=0, pz_4jets=0, pt_norm=0, DR_HiggsVBF=0, Dphi_qq_cc=0, jetEne_sum=0, jetPt_sum=0;
         //   int njets=0;
         //
         // CvsAll_jetC1, CvsAll_jetC2, CvsB_jetC1, CvsB_jetC2, CvsL_jetC1, CvsL_jetC2
         CvsAll_jetC1=AK4puppi_cHiggs.at(0).CvsAll;
         CvsAll_jetC2=AK4puppi_cHiggs.at(1).CvsAll;
         
         CvsL_jetC1=AK4puppi_cHiggs.at(0).CvsL;
         CvsL_jetC2=AK4puppi_cHiggs.at(1).CvsL;
         CvsB_jetC1=AK4puppi_cHiggs.at(0).CvsB;
         CvsB_jetC2=AK4puppi_cHiggs.at(1).CvsB;
         mqq = (jet_VBF1+jet_VBF2).M();
         Deta_qq = fabs(jet_VBF1.Eta()-jet_VBF2.Eta());
         Dphi_qq = fabs(jet_VBF1.DeltaPhi(jet_VBF2));
         TLorentzVector VBFsum= jet_VBF1+jet_VBF2;
         TVector3 VBFsum_vec;
         VBFsum_vec.SetPtEtaPhi(VBFsum.Pt(), VBFsum.Eta(), VBFsum.Phi());
         float alpha_VBF1 = jet_VBF1.Angle(VBFsum_vec);
         float alpha_VBF2 = jet_VBF2.Angle(VBFsum_vec);
         Alpha_qq = std::min(alpha_VBF1, alpha_VBF2);
         qgl_VBF1 = AK4puppi_VBF.at(0).QGL;
         qgl_VBF2 = AK4puppi_VBF.at(1).QGL;
         pz_4jets = jet_c1.Pz() + jet_c2.Pz() + jet_VBF1.Pz() + jet_VBF2.Pz();
         TLorentzVector p4_4jets = jet_c1 + jet_c2 + jet_VBF1 + jet_VBF2;
         pt_norm = p4_4jets.Pt()/(jet_c1.Pt() + jet_c2.Pt() + jet_VBF1.Pt() + jet_VBF2.Pt());
         TLorentzVector Higgs_p4 = jet_c1 + jet_c2;
         DR_HiggsVBF1 = Higgs_p4.DeltaR(jet_VBF1);
         DR_HiggsVBF2 = Higgs_p4.DeltaR(jet_VBF2);
         Dphi_qq_cc = fabs(VBFsum.DeltaPhi(Higgs_p4));
         mCC= (jet_c1+jet_c2).M(); 
         tree->Fill();
       }
     }//end condition on HLT
     for(int i=0; i<NCUTS; i++){
       if(Eff_counter[i] == true) cutevt[i]++;
     }
 } //end cycle on events



    TCanvas *canvEvt = new TCanvas("CutEfficiency_Nevents", "CutEfficiency_Nevents", 0, 0, 1200, 1000);
    Draw_CutEffCanvas(canvEvt, hCutEffEvt, cutevt, listCut);
    
    //Write and close the file
   
   histo_pt_jet1->Write();
   histo_pt_jet2->Write();
   histo_pt_jet3->Write();
   histo_pt_jet4->Write();
   histo_eta_jet1->Write();
   histo_eta_jet2->Write();
   histo_eta_jet3->Write();
   histo_eta_jet4->Write();
   histo_pt_jetC1->Write();
   histo_pt_jetC2->Write();
   histo_pt_jetVBF1->Write();
   histo_pt_jetVBF2->Write();
   histo_eta_jetC1->Write();
   histo_eta_jetC2->Write();
   histo_eta_jetVBF1->Write();
   histo_eta_jetVBF2->Write();
   histo_deltaEta_cjets->Write();
   histo_deltaEta_VBFjets->Write();
   histo_invMass_cjets->Write();
   histo_invMass_VBFjets->Write();
   histo_CvsAll_c1->Write();
   histo_CvsAll_c2->Write();
   histo_CvsL_c1->Write();
   histo_CvsL_c2->Write();
   histo_CvsB_c1->Write();
   histo_CvsB_c2->Write();

    fout->Write();
    fout->Close();
    

    /*TFile *rootFile = new TFile("AK4puppi_selEvents.root","RECREATE");
    CvsAll_ptSel->Write();
    pt_cjets->Write();
    pt_cjets->Write();
    pt_VBFjets->Write();
    pt_VBFjets->Write();
    eta_cjets->Write();
    eta_cjets->Write();
    eta_VBFjets->Write();
    eta_VBFjets->Write();
    deltaEta_cjets->Write();
    deltaEta_VBFjets->Write();
    invMass_cjets->Write();
    invMass_VBFjets->Write();
    CvsAll_c1->Write();
    CvsAll_c2->Write();
    CvsL_c1->Write();
    CvsL_c2->Write();
    CvsB_c1->Write();
    CvsB_c2->Write();
    rootFile->Close();*/
}
