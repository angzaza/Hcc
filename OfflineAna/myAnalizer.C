#define myAnalizer_cxx
#define NCUTS 10
#include "myAnalizer.h"
#include "Utilities.C"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>


using namespace std;

void Jet::print(){
  cout<<"pt: "<<pt_<<"  eta: "<<eta_<<"  phi: "<<phi_<<"   mass: "<<mass_<<endl;
  cout<<"Cvsll: "<<CvsAll_<<"   CvsL"<<CvsL_<<"   CvsB: "<<CvsB_<<endl;
  cout<<"QvsG: "<<QvsG_<<"   QGL: "<<QGL_<<endl; 
  cout<<"ctagged: "<<ctagged_<<endl;
}

bool compareByPt(const TLorentzVector &a, const TLorentzVector &b){
  return a.Pt() > b.Pt();
}

/*struct jet_struct{
  TLorentzVector jet_p4;
  float CvsAll;
  float CvsL;
  float CvsB;
  float QGL;
  float QvsG;
  bool c_tagged;
};*/

bool compareByCvsAll(const Jet &a, const Jet &b){
  if (std::isnan(a.CvsAll())) {
    // If a is NaN, b should come before a
    return false;
  } 
  else if (std::isnan(b.CvsAll())) {
    // If b is NaN, a should come before b
    return true;
  } 
  else{
     // Neither a nor b is NaN, compare their ctag values
     return a.CvsAll() > b.CvsAll();
  }
}

bool compareByCvsL(const Jet &a, const Jet &b){
  if (std::isnan(a.CvsL())) {
    // If a is NaN, b should come before a
    return false;
  } 
  else if (std::isnan(b.CvsL())) {
        // If b is NaN, a should come before b
        return true;
  } 
  else {
    // Neither a nor b is NaN, compare their ctag values
    return a.CvsL() > b.CvsL();
  }
}

bool IsVetoed(TH2D* histo, double eta, double phi){
  int binX = histo->GetXaxis()->FindBin(eta);
  int binY = histo->GetYaxis()->FindBin(phi);
  if (histo->GetBinContent(binX, binY) != 0) return true;
  else return false;
}

float Get_PUweight(TH1F* histo, int pu){
  int binX = histo->GetXaxis()->FindBin(pu);
  float puw = histo->GetBinContent(binX);
  //cout<<"apply veto function cc"<<endl;
    return puw;
}


void myAnalizer::Loop_Hcc(TString type, TString datasetName, TString era)
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

   bool verbose=false;

   //variable and histo declaration
   TLorentzVector hltjet_p4;
   TLorentzVector GENjet_p4;
   TLorentzVector AK4puppi_p4;
   //vector<TLorentzVector> AK4puppi_p4_vec;
   //vector<TLorentzVector> cjets_hlt, VBFjets_hlt, fakejets_hlt;
   TVector3 hltjet_p3, quark_p3, GENjet_p3;
   int match_index=-10;
   int match_hltGen_index=-10;
   vector<bool> match_free, match_gen_free;
   bool match_gen_quark=false;
   bool match_hlt_gen=false;
   bool trigger_VBFPNet=false;
   bool Apply_jetVeto=true;
   bool Apply_PUrew=true;
   
   vector<Jet> goodAK4puppiJet;
   vector<float> hltjet_pt;

   int hlt_count=0;
   int hltdouble_count=0;
   int hltOR_count=0;

   bool PtThr_passed=false;
   float pt_hltThr[3]={100.,88.,70.};
 
   //define a vector of jet_struct
   //std::vector<jet_struct> AK4puppi_sel;
   //std::vector<jet_struct> AK4puppi_cHiggs;
   //std::vector<jet_struct> AK4puppi_VBF;
   std::vector<Jet> AK4puppi_sel;
   std::vector<Jet> AK4puppi_cHiggs;
   std::vector<Jet> AK4puppi_VBF;
   bool passed_selection=false;
   //
   //
   //
   // integrated lumi in fb ////////////
   float Lumi_2023C=17.981;
   

   // Open CSV file for writing
   std::ofstream outputFile("VBFHCC.csv");
   outputFile <<"event,jet1pt,jet1eta,jet1phi,jet1ene,jet1CvsAll,jet2pt,jet2eta,jet2phi,jet2ene,jet2CvsAll,jet3pt,jet3eta,jet3phi,jet3ene,jet3CvsAll,jet4pt,jet4eta,jet4phi,jet4ene,jet4CvsAll,"<<std::endl;

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
   TH1F *histo_pt_jetC1=new TH1F("pt_jetc1","pt_jetc1",25,30,330);
   TH1F *histo_pt_jetC2=new TH1F("pt_jetc2","pt_jetc2",25,30,330);
   TH1F *histo_pt_jetVBF1=new TH1F("pt_jetVBF1","pt_jetVBF1",25,30,330);
   TH1F *histo_pt_jetVBF2=new TH1F("pt_jetVBF2","pt_jetVBF2",25,30,330);
   TH1F *histo_eta_jetC1=new TH1F("eta_jetc1"," eta_jetc1",50,-5,5);
   TH1F *histo_eta_jetC2=new TH1F("eta_jetc2"," eta_jetc2",50,-5,5);
   TH1F *histo_eta_jetVBF1=new TH1F("eta_jetVBF1"," eta_jetVBF2",50,-5,5);
   TH1F *histo_eta_jetVBF2=new TH1F("eta_jetVBF2"," eta_jetVBF2",50,-5,5);
   histo_pt_jet1->Sumw2();
   histo_pt_jet2->Sumw2();
   histo_pt_jet3->Sumw2();
   histo_pt_jet4->Sumw2();
   histo_eta_jet1->Sumw2();
   histo_eta_jet2->Sumw2();
   histo_eta_jet3->Sumw2();
   histo_eta_jet4->Sumw2();
   histo_pt_jetC1->Sumw2();
   histo_pt_jetC2->Sumw2();
   histo_pt_jetVBF1->Sumw2();
   histo_pt_jetVBF2->Sumw2();
   histo_eta_jetC1->Sumw2();
   histo_eta_jetC2->Sumw2();
   histo_eta_jetVBF1->Sumw2();
   histo_eta_jetVBF2->Sumw2();

   TH1F* histo_deltaEta_cjets= new TH1F("DeltaEta_cjets","DelatEta_cjets",50,0,12);
   TH1F* histo_deltaEta_VBFjets= new TH1F("DeltaEta_VBFjets","DeltaEta_VBFjets",50,2.2,7.2);
   histo_deltaEta_cjets->Sumw2();
   histo_deltaEta_VBFjets->Sumw2();

   TH1F* histo_invMass_cjets= new TH1F("InvMass_cjets","InvMass_cjets",50,50,200);
   TH1F* histo_invMass_VBFjets= new TH1F("InvMass_VBFjets","InvMass_VBFjets",50,400,3000);
   histo_invMass_cjets->Sumw2();
   histo_invMass_VBFjets->Sumw2();

   TH1F* histo_CvsAll_c1= new TH1F("PNet_CvsAll_jetc1","PNet_CvsAll_jetc1",50,0,1);
   TH1F* histo_CvsAll_c2= new TH1F("PNet_CvsAll_jetc2","PNet_CvsAll_jetc2",50,0,1);
   TH1F* histo_CvsL_c1= new TH1F("PNet_CvsL_jetc1","PNet_CvsL_jetc1",50,0,1);
   TH1F* histo_CvsL_c2= new TH1F("PNet_CvsL_jetc2","PNet_CvsL_jetc2",50,0,1);
   TH1F* histo_CvsB_c1= new TH1F("PNet_CvsB_jetc1","PNet_CvsB_jetc1",50,0,1);
   TH1F* histo_CvsB_c2= new TH1F("PNet_CvsB_jetc2","PNet_CvsB_jetc2",50,0,1);
   histo_CvsAll_c1->Sumw2();
   histo_CvsAll_c2->Sumw2();
   histo_CvsL_c1->Sumw2();
   histo_CvsL_c2->Sumw2();
   histo_CvsB_c1->Sumw2();
   histo_CvsB_c1->Sumw2();

   TH1F* histo_QvsG_VBF1= new TH1F("PNet_QvsG_jetVBF1","PNet_QvsG_jetVBF1",50,-1,1);
   TH1F* histo_QvsG_VBF2= new TH1F("PNet_QvsG_jetVBF2","PNet_QvsG_jetVBF2",50,1,1);
   histo_QvsG_VBF1->Sumw2();
   histo_QvsG_VBF1->Sumw2();

   double isMC = -99;
   double run_n = 0, lumi_n = 0, evt_n = 0, pileupFactor=0;
   double XS = 1.;
   double XS_err = 1.;
   bool HLT_passed=false;
   bool HLT_VBFHCC=false;
   bool HLT_parking=false;
   double Lumi=-1.;

   double pt_jetC1=0, pt_jetC2=0, pt_jetVBF1=0, pt_jetVBF2=0, eta_jetC1=0, eta_jetC2=0, eta_jetVBF1=0, eta_jetVBF2=0, CvsAll_jetC1=0, CvsAll_jetC2=0, CvsB_jetC1=0, CvsB_jetC2=0, CvsL_jetC1=0, CvsL_jetC2=0;
   double mqq=0, Deta_qq=0, Dphi_qq=0, Alpha_qq=0, qgl_VBF1=0, qgl_VBF2=0, QvsG_VBF1=0, QvsG_VBF2=0, pz_4jets=0, pt_norm=0, DR_HiggsVBF1=0, DR_HiggsVBF2=0,  Dphi_qq_cc=0, jetEne_sum=0, jetPt_sum=0, mCC=0;
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
   TreeFin_Init(tree, isMC, lumi_n, run_n, evt_n, entry, XS, XS_err, Lumi, pileupFactor, pt_jetC1, pt_jetC2, pt_jetVBF1, pt_jetVBF2, eta_jetC1, eta_jetC2, eta_jetVBF1, eta_jetVBF2, CvsAll_jetC1, CvsAll_jetC2, CvsB_jetC1, CvsB_jetC2, CvsL_jetC1, CvsL_jetC2, mqq, Deta_qq, Dphi_qq, Alpha_qq, qgl_VBF1, qgl_VBF2, QvsG_VBF1, QvsG_VBF2, pz_4jets, pt_norm, DR_HiggsVBF1, DR_HiggsVBF2  ,Dphi_qq_cc, njets, jetEne_sum, jetPt_sum, mCC);

   //name of the jet veto map file
   TString vetoMapFile;
   TFile* vetoFile;
   TH2D* vetoMapHisto;

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


  // PU reweighting
   TFile* PUrewFile;
   TH1F* PUrewHisto;
 
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

   

   if(datasetName.Contains("2023")) isMC=0;
   else{
        if(datasetName.Contains("VBFHCC")) isMC=1;
        if(datasetName.Contains("ggHcc")) isMC=2;
        if(datasetName.Contains("VBFHBB")) isMC=3;
        if(datasetName.Contains("ggHBB")) isMC=4;
    }

   // initialize Cross Section
   if(isMC!=0){
     XS = CrossSection(datasetName).first;
     XS_err = CrossSection(datasetName).second;
   }

   //if(datasetName.Contains("2023C")) Lumi=Lumi_2023C; //fb
   if(era.Contains("2023C")) Lumi=Lumi_2023C; //fb

   if (fChain == 0) return;

   //Long64_t nentries = 1000000;
   Long64_t nentries = fChain->GetEntriesFast();
   //Long64_t nentries = 1000;

   Long64_t nbytes = 0, nb = 0;
   
   //cycle on events
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      entry=jentry;
      if(verbose) cout<<"evt: "<<entry<<endl;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      bool Eff_counter[NCUTS] = {false};
      goodAK4puppiJet.clear();
      AK4puppi_sel.clear();
      AK4puppi_cHiggs.clear();
      AK4puppi_VBF.clear();
      HLT_passed=false;
      HLT_VBFHCC=false;
      HLT_parking=false;
    
      //set default values
      PtThr_passed=false;
      passed_selection=false;

      int n_hltmatch=0;

      Eff_counter[0] = true;

      run_n = Run; lumi_n = LumiSect; evt_n = Event; pileupFactor=nInt;  
    
      njets=0; jetEne_sum=0; jetPt_sum=0;
 
      //HLT trigger
      for(int h=0; h<Trigger_hltname->size(); h++){
        TString hltName = Trigger_hltname->at(h);

         if( (hltName.Contains("HLT_QuadPFJet100_88_70_30_PNet1CvsAll0p5_VBF3Tight_v")) && Trigger_hltdecision->at(h) == 1){
          HLT_passed = true;
          HLT_VBFHCC=true;
          //Eff_counter[1] = true;
         }


         //if( (hltName.Contains("HLT_VBF_DiPFJet105_40_Mjj1000_Detajj3p5_v")) && Trigger_hltdecision->at(h) == 1){
         //  HLT_parking=true;
         //  HLT_passed=true;
         //}
      }

      if(!HLT_passed){ // se non passa HLT, aggiorna il conteggio degli eventi e passa all'evento successivo
        cutevt[0]++;
        continue;
      }

      //Pileup reweighting
      float pu_weight=1.0;
      if(isMC!=0 && Apply_PUrew==true) pu_weight=Get_PUweight(PUrewHisto, nInt);

      //if(verbose){for(unsigned int i=0; i<quark_VBF->size(); i++){ cout<<"quark VBF: "<<quark_VBF->at(i)<<endl;}}

      if(HLT_passed==true){
      //if(HLT_VBFHCC==true && HLT_parking==false){
        Eff_counter[1] = true;

        bool LeptonVeto=false;
        //veto on electrons and muons
        for(unsigned int iele=0; iele< Ele_pt->size(); iele++){
          if(Ele_pt->at(iele)<7.0 || fabs(Ele_eta->at(iele)>2.5) || fabs(Ele_dxy->at(iele))>0.05 || fabs(Ele_dz->at(iele))>0.2) continue;
          if(Ele_IsoCal->at(iele)<=0.4){
            LeptonVeto=true;
            break;
          }
        }
        if(LeptonVeto==false){
          for(unsigned int imu=0; imu< Muon_pt->size(); imu++){
            if(Muon_pt->at(imu)<5.0 || fabs(Muon_eta->at(imu)>2.5) || fabs(Muon_dxy->at(imu))>0.25 || fabs(Muon_dz->at(imu))>1.0) continue;
            if(Muon_PF_Iso_R04->at(imu)<=0.4){
              LeptonVeto=true;
              break;
            }
          }
        } 
        if(LeptonVeto==false && met_pt<170.) Eff_counter[2] = true;
        if(LeptonVeto==true || met_pt>170.){
          for(int i=0; i<NCUTS; i++){
            if(Eff_counter[i] == true) cutevt[i]++;
          }
          continue; //skip current event and update the counters 

        }
          //cout<<"veto false"<<endl;
          //Eff_counter[2] = true;

         bool JetVetoed=false;

          // check offline pt threshold first 4 jets
          if(AK4PuppiJets_pt->size() < 4){
             for(int i=0; i<NCUTS; i++){
               if(Eff_counter[i] == true) cutevt[i]++;
             }
             continue; //skip current event and update counter
          }
          float pt0= AK4PuppiJets_pt->at(0);
          float pt1= AK4PuppiJets_pt->at(1);
          float pt2= AK4PuppiJets_pt->at(2);
          float pt3= AK4PuppiJets_pt->at(3);
          if(verbose) cout<<"pt: "<<pt0<<"  "<<pt1<<"  "<<pt2<<"  "<<pt3<<endl;
          if(pt0<105. || pt1<90. || pt2<75. || pt3<35.){
            for(int i=0; i<NCUTS; i++){
              if(Eff_counter[i] == true) cutevt[i]++;
            }
            continue;
          }
          Eff_counter[3] = true; // counter offline pt theshold
          if(verbose) cout<<"pt offline threshold passed "<<endl;
          //cycle on AK4 puppi jets
          for(unsigned int ijet=0; ijet< AK4PuppiJets_pt->size();ijet++){
          //for(unsigned int ijet=0; ijet<4;ijet++){
            TVector3 AK4puppi_p3;

            float num= jet_pfParticleNetAK4JetTags_probc->at(ijet);
            float den= jet_pfParticleNetAK4JetTags_probc->at(ijet)+ jet_pfParticleNetAK4JetTags_probb->at(ijet) + jet_pfParticleNetAK4JetTags_probuds->at(ijet) + jet_pfParticleNetAK4JetTags_probg->at(ijet);
            //float CvsAll_val= (den!=0 || (num/den)<999)? num/den : -1;
            float CvsAll_val=num/den;
            float QGL_val=AK4PuppiJets_qgl->at(ijet);
            float CvsB_val=jet_pfParticleNetAK4JetTags_CvsB->at(ijet);
            float CvsL_calc = jet_pfParticleNetAK4JetTags_probc->at(ijet)/(jet_pfParticleNetAK4JetTags_probc->at(ijet)+jet_pfParticleNetAK4JetTags_probuds->at(ijet) + jet_pfParticleNetAK4JetTags_probg->at(ijet));
            float CvsL_val=jet_pfParticleNetAK4JetTags_CvsL->at(ijet);
            //cout<<"CvsL calc: "<<CvsL_calc<<"   CvsL def: "<<CvsL_val<<endl;
            float QvsG_val=jet_pfParticleNetAK4JetTags_QvsG->at(ijet);
            Jet jet_i(AK4PuppiJets_pt->at(ijet),AK4PuppiJets_eta->at(ijet), AK4PuppiJets_phi->at(ijet), AK4PuppiJets_mass->at(ijet), CvsAll_val, CvsB_val, CvsL_val, QvsG_val, QGL_val, false);
            if(verbose) cout<<"#############################################################"<<endl;
            if(verbose) cout<<"jet "<<ijet<<endl;
            if(verbose) jet_i.print();
            AK4puppi_p3.SetPtEtaPhi(AK4PuppiJets_pt->at(ijet), AK4PuppiJets_eta->at(ijet), AK4PuppiJets_phi->at(ijet));
            //check if the jet contains a muon in dR<0.2
            /*float dRjmu_max=0.2;
            bool match_jmu=false;
            if(verbose) cout<<"matching with muons"<<endl;
            for(unsigned int j=0;j<Muon_pt->size();j++){
               TVector3 muon_p3;
               muon_p3.SetPtEtaPhi(Muon_pt->at(j),Muon_eta->at(j), Muon_phi->at(j));
               float dRjmu = AK4puppi_p3.DeltaR(muon_p3); 
               if(verbose) cout<<"dR:"<< dRjmu<<endl;
               if(dRjmu<=dRjmu_max){
                 match_jmu=true;
                 if(verbose) cout<<"matched -- skipping current jet"<<endl;
                 break;
               } 
            }
            //if(match_jmu==true) continue; 
            if(match_jmu==true && ijet<4) break; //se uno dei primi 4 jet contiene un muone devo scartare l'evento
            else if(match_jmu==true && ijet>=4) continue; //se un altro dei jet dell'evento proviene da un muone, scarto solo il jet 
           */
            //goodAK4puppiJet.push_back(jet);

            //if a jet is vetoed by the jet veto map, the event is discarded 
            if(Apply_jetVeto==true){
              //for(unsigned int ijet=0, ijet<goodAK4puppiJet.size();ijet++){
              if(IsVetoed(vetoMapHisto,AK4PuppiJets_eta->at(ijet), AK4PuppiJets_phi->at(ijet) )){
                //if(!ApplyVeto(vetoMapHisto,goodAK4puppiJet[ijet].eta(), goodAK4puppiJet[ijet].phi() )){
                JetVetoed=true;
                if(verbose) cout<<"event vetoed by jet maps  -- skipping current event"<<endl;
                break;
              }
            }
            //if(ijet==AK4PuppiJets_eta->size()-1) Eff_counter[4]=true;
            //if(JetVetoed==false) Eff_counter[2]=true;
            bool match_hltjet=false;
            if(ijet<3){  
              //TVector3 AK4puppi_p3;
              //AK4puppi_p3.SetPtEtaPhi(AK4PuppiJets_pt->at(ijet), AK4PuppiJets_eta->at(ijet), AK4PuppiJets_phi->at(ijet));
              //match with L3 trigger objects
              if(verbose) cout<<"match with HLT jets: "<<endl;
              float dR_max=0.4;
              for(unsigned int j=0; j<HLTAK4PFJetLoose_pt->size(); j++){
                TVector3 hltjet_p3;
                hltjet_p3.SetPtEtaPhi(HLTAK4PFJetLoose_pt->at(j), HLTAK4PFJetLoose_eta->at(j), HLTAK4PFJetLoose_phi->at(j));
                if(verbose) cout<<"hlt jet "<<j<<"  pt: "<<HLTAK4PFJetLoose_pt->at(j)<<endl;
                if( hltjet_p3.Pt() < pt_hltThr[j] ) continue;
                float this_dR=AK4puppi_p3.DeltaR(hltjet_p3);
                if(verbose) cout<<"dR:"<<this_dR<<endl;
                if(this_dR<dR_max){
                  match_hltjet=true;
                  n_hltmatch++;
                  if(verbose) cout<<"matched -- breaking current cycle"<<endl;
                  break;
                }
              }
             }       
            
            if(match_hltjet==false && ijet<3){
               break; //se uno dei primi 3 jet non matcha con un hlt object esco dal ciclo dei jet
               if(verbose) cout<<" one of the first three offline jets doesn't match with an hlt object -- skipping current event"<< endl;
            }
            //if(match_hltjet==true){
            if(fabs(AK4puppi_p4.Eta())<4.7) AK4puppi_sel.push_back(jet_i); 

            if(AK4puppi_p4.Pt()>20. && fabs(AK4puppi_p4.Eta()<2.4)) njets++;
          }//end cycle on puppi
          if(JetVetoed==true){
            for(int i=0; i<NCUTS; i++){
              if(Eff_counter[i] == true) cutevt[i]++;
            }
            continue;
          }
          Eff_counter[4] = true;
          if(n_hltmatch!=3){
            for(int i=0; i<NCUTS; i++){
              if(Eff_counter[i] == true) cutevt[i]++;
            }
            continue;
          }
          Eff_counter[5] = true;
         //if(JetVetoed==false){
          //  Eff_counter[2]=true;
          if(verbose) cout<<"##############################"<<endl;
          if(verbose) cout<<"size AK4puppi_sel "<<AK4puppi_sel.size()<<endl;
          if(AK4puppi_sel.size()>=4){
            for(unsigned int j=4; j<AK4puppi_sel.size(); j++){
              TLorentzVector AK4puppi_sel_p4;
              AK4puppi_sel_p4.SetPtEtaPhiM(AK4puppi_sel[j].pt(),AK4puppi_sel[j].eta(), AK4puppi_sel[j].phi(), AK4puppi_sel[j].mass());
              if(AK4puppi_sel_p4.Pt()<30 || fabs(AK4puppi_sel_p4.Eta())>2.4) continue; 
              jetEne_sum+=AK4puppi_sel_p4.E();
              jetPt_sum+=AK4puppi_sel_p4.Pt(); 
          }
          //resize the vector of selected jets, keeping only the first 4 of them (pt sorted)
          AK4puppi_sel.resize(4);
          //cout<<"AK4puppi_sel resized at 4"<<endl;
          float pt0= AK4puppi_sel[0].pt();
          float pt1= AK4puppi_sel[1].pt();
          float pt2= AK4puppi_sel[2].pt();
          float pt3= AK4puppi_sel[3].pt();
          if(verbose) cout<<"pt: "<<pt0<<"  "<<pt1<<"  "<<pt2<<"  "<<pt3<<endl;
          //if(pt0>=105.0 && pt1>= 90. && pt2>=75. && pt3>=35.){
          //  if(verbose) cout<<"pt threshold passed"<<endl;
          //PtThr_passed=true;
          histo_pt_jet1->Fill(pt0, pu_weight);
          histo_pt_jet2->Fill(pt1, pu_weight);
          histo_pt_jet3->Fill(pt2, pu_weight);
          histo_pt_jet4->Fill(pt3, pu_weight);
          histo_eta_jet1->Fill(AK4puppi_sel[0].eta(), pu_weight);
          histo_eta_jet2->Fill(AK4puppi_sel[1].eta(), pu_weight);
          histo_eta_jet3->Fill(AK4puppi_sel[2].eta(), pu_weight);
          histo_eta_jet4->Fill(AK4puppi_sel[3].eta(), pu_weight);
          //}
          //if(PtThr_passed==true){
          //Eff_counter[3]=true;
          //sort the AK4puppi_sel by descending CvsAll score
          std::sort(AK4puppi_sel.begin(),AK4puppi_sel.end(),compareByCvsAll );
          // CvsAll_ptSel->Fill(AK4puppi_sel.at(0).CvsAll);
          //tag as c jets the two jets with the highest CvsAll score within eta<2.4
          int n_cjets=0;
          if(verbose) cout<<" checking if 2 jets are tagged as c: "<<endl;
            for(unsigned int k=0; k<4; k++){
					   if(n_cjets<2){ //DA CONTROLLARE
                if(fabs((AK4puppi_sel[k].eta()))<=2.4){
                AK4puppi_sel[k].set_ctagged(true);
                if(verbose) cout<<"jet "<<k<<"  eta: "<<AK4puppi_sel[k].eta()<<"  ctagged: "<<AK4puppi_sel[k].ctagged()<<endl;
                n_cjets++;
              }
            }  
          }
          if(n_cjets==2){
            Eff_counter[6]=true;
            if(verbose) cout<<"2 c-jets found"<<endl;
            for(unsigned int j=0; j<4; j++){
              if(AK4puppi_sel[j].ctagged()==true){
                AK4puppi_cHiggs.push_back(AK4puppi_sel[j]);
                if(verbose) cout<<"c jet:"<<endl;
                if(verbose) AK4puppi_sel[j].print();
              }
              else{
                AK4puppi_VBF.push_back(AK4puppi_sel[j]);
                if(verbose) cout<<"VBF jet"<<endl;
                if(verbose) AK4puppi_sel[j].print();
              }
            }
            if(AK4puppi_cHiggs.size()==2 && AK4puppi_VBF.size()==2){
              Eff_counter[7]=true;
              if(verbose){if(n_cjets!=2){cout<<"AAAAAAAA"<<endl;}} 
              if(AK4puppi_cHiggs[0].CvsAll()>=0.5 && AK4puppi_cHiggs[0].CvsL()>=0.16 && AK4puppi_cHiggs[0].CvsB()>=0.304 
                 && AK4puppi_cHiggs[1].CvsL()>=0.054 && AK4puppi_cHiggs[1].CvsB()>=0.182 ){ //condition on ctag verified -- Medium WPs for 2022E (to be updated)
                Eff_counter[8]=true;
                if(verbose) cout<<"passed c tag selection"<<endl;
                TLorentzVector jet_VBF1;
                jet_VBF1.SetPtEtaPhiM(AK4puppi_VBF[0].pt(),AK4puppi_VBF[0].eta(), AK4puppi_VBF[0].phi(), AK4puppi_VBF[0].mass());
                TLorentzVector jet_VBF2;
                jet_VBF2.SetPtEtaPhiM(AK4puppi_VBF[1].pt(),AK4puppi_VBF[1].eta(), AK4puppi_VBF[1].phi(), AK4puppi_VBF[1].mass());
                if(verbose) cout<<"mass VBF jets: "<<(jet_VBF1+jet_VBF2).M()<<"    Deta: "<<fabs(jet_VBF1.Eta()-jet_VBF2.Eta())<<endl;
                TLorentzVector jet_cHiggs1;
                jet_cHiggs1.SetPtEtaPhiM(AK4puppi_cHiggs[0].pt(),AK4puppi_cHiggs[0].eta(), AK4puppi_cHiggs[0].phi(), AK4puppi_cHiggs[0].mass());
                TLorentzVector jet_cHiggs2;
                jet_cHiggs2.SetPtEtaPhiM(AK4puppi_cHiggs[1].pt(),AK4puppi_cHiggs[1].eta(), AK4puppi_cHiggs[1].phi(), AK4puppi_cHiggs[1].mass());
                if(verbose) cout<<"mass c jets: "<<(jet_cHiggs1+jet_cHiggs2).M()<<"    Deta: "<<fabs(jet_cHiggs1.Eta()-jet_cHiggs2.Eta())<<endl;
                if((jet_VBF1+jet_VBF2).M()>=500. && fabs(jet_VBF1.Eta()-jet_VBF2.Eta())>=3.8){
                  passed_selection=true;
                   if(verbose) cout<<"passed selection"<<endl;
                }
              }
            }
          } 
          //}
        }
        // }
       //}

       //PassedEvents->Fill(passed_selection);

       if(passed_selection==true){
         Eff_counter[9]=true;
         //TLorentzVector jet_VBF1=AK4puppi_VBF.at(0).jet_p4;
         //TLorentzVector jet_VBF2=AK4puppi_VBF.at(1).jet_p4;
         //TLorentzVector jet_c1=AK4puppi_cHiggs.at(0).jet_p4;
         //TLorentzVector jet_c2=AK4puppi_cHiggs.at(1).jet_p4;

         TLorentzVector jetVBF1_p4;
         jetVBF1_p4.SetPtEtaPhiM(AK4puppi_VBF[0].pt(),AK4puppi_VBF[0].eta(), AK4puppi_VBF[0].phi(), AK4puppi_VBF[0].mass());
         TLorentzVector jetVBF2_p4;
         jetVBF2_p4.SetPtEtaPhiM(AK4puppi_VBF[1].pt(),AK4puppi_VBF[1].eta(), AK4puppi_VBF[1].phi(), AK4puppi_VBF[1].mass());
         TLorentzVector jetc1_p4;
         jetc1_p4.SetPtEtaPhiM(AK4puppi_cHiggs[0].pt(),AK4puppi_cHiggs[0].eta(), AK4puppi_cHiggs[0].phi(), AK4puppi_cHiggs[0].mass());
         TLorentzVector jetc2_p4;
         jetc2_p4.SetPtEtaPhiM(AK4puppi_cHiggs[1].pt(),AK4puppi_cHiggs[1].eta(), AK4puppi_cHiggs[1].phi(), AK4puppi_cHiggs[1].mass());

         histo_pt_jetC1->Fill(AK4puppi_cHiggs[0].pt(), pu_weight);
         histo_pt_jetC2->Fill(AK4puppi_cHiggs[1].pt(), pu_weight);

         histo_eta_jetC1->Fill(AK4puppi_cHiggs[0].eta(), pu_weight);
         histo_eta_jetC2->Fill(AK4puppi_cHiggs[1].eta(), pu_weight);

         histo_pt_jetVBF1->Fill(AK4puppi_VBF[0].pt(), pu_weight);
         histo_pt_jetVBF2->Fill(AK4puppi_VBF[1].pt(),pu_weight);

         histo_eta_jetVBF1->Fill(AK4puppi_VBF[0].eta(),pu_weight);
         histo_eta_jetVBF2->Fill(AK4puppi_VBF[1].eta(),pu_weight);
         histo_QvsG_VBF1->Fill(AK4puppi_VBF[0].QvsG()),pu_weight;
         histo_QvsG_VBF2->Fill(AK4puppi_VBF[1].QvsG()),pu_weight;
         histo_deltaEta_cjets->Fill(fabs(AK4puppi_cHiggs[0].eta()-AK4puppi_cHiggs[1].eta()),pu_weight);
         histo_deltaEta_VBFjets->Fill(fabs(AK4puppi_VBF[1].eta()-AK4puppi_VBF[1].eta()),pu_weight);

         histo_invMass_cjets->Fill((jetc1_p4+jetc2_p4).M(),pu_weight);
         histo_invMass_VBFjets->Fill((jetVBF1_p4+jetVBF2_p4).M(),pu_weight);

         histo_CvsAll_c1->Fill(AK4puppi_cHiggs[0].CvsAll(),pu_weight);
         histo_CvsAll_c2->Fill(AK4puppi_cHiggs[1].CvsAll(),pu_weight);

         histo_CvsL_c1->Fill(AK4puppi_cHiggs[0].CvsL(),pu_weight);
         histo_CvsL_c2->Fill(AK4puppi_cHiggs[1].CvsL(),pu_weight);

         histo_CvsB_c1->Fill(AK4puppi_cHiggs[0].CvsB(),pu_weight);
         histo_CvsB_c2->Fill(AK4puppi_cHiggs[1].CvsB(),pu_weight);

         pt_jetC1 = AK4puppi_cHiggs[0].pt();
         pt_jetC2 = AK4puppi_cHiggs[1].pt();
          
         eta_jetC1 = AK4puppi_cHiggs[0].eta();
         eta_jetC2 = AK4puppi_cHiggs[1].eta();

         pt_jetVBF1 = AK4puppi_VBF[0].pt();
         pt_jetVBF2 = AK4puppi_VBF[1].pt();

         eta_jetVBF1 = AK4puppi_VBF[0].eta();
         eta_jetVBF2 = AK4puppi_VBF[1].eta();

         // BDT variables
         //Deta_qq=0, Dphi_qq=0, Alfa_qq=0, qgl_qq=0, pz_4jets=0, pt_norm=0, DR_HiggsVBF=0, Dphi_qq_cc=0, jetEne_sum=0, jetPt_sum=0;
         //   int njets=0;
         //
         // CvsAll_jetC1, CvsAll_jetC2, CvsB_jetC1, CvsB_jetC2, CvsL_jetC1, CvsL_jetC2
         CvsAll_jetC1=AK4puppi_cHiggs[0].CvsAll();
         CvsAll_jetC2=AK4puppi_cHiggs[1].CvsAll();
         
         CvsL_jetC1=AK4puppi_cHiggs[0].CvsL();
         CvsL_jetC2=AK4puppi_cHiggs[1].CvsL();
         CvsB_jetC1=AK4puppi_cHiggs[0].CvsB();
         CvsB_jetC2=AK4puppi_cHiggs[1].CvsB();
         mqq = (jetVBF1_p4+jetVBF2_p4).M();
         Deta_qq = fabs(fabs(jetVBF1_p4.Eta()-jetVBF2_p4.Eta()));
         Dphi_qq = fabs(jetVBF1_p4.DeltaPhi(jetVBF2_p4));
         TLorentzVector VBFsum= jetVBF1_p4+jetVBF2_p4;
         TVector3 VBFsum_vec;
         VBFsum_vec.SetPtEtaPhi(VBFsum.Pt(), VBFsum.Eta(), VBFsum.Phi());
         TVector3 jetVBF1_p3;
         jetVBF1_p3.SetPtEtaPhi(jetVBF1_p4.Pt(), jetVBF1_p4.Eta(), jetVBF1_p4.Phi());
         TVector3 jetVBF2_p3;
         jetVBF2_p3.SetPtEtaPhi(jetVBF2_p4.Pt(), jetVBF2_p4.Eta(), jetVBF2_p4.Phi());
         float alpha_VBF1 = jetVBF1_p3.Angle(VBFsum_vec);
         float alpha_VBF2 = jetVBF2_p3.Angle(VBFsum_vec);
         Alpha_qq = std::min(alpha_VBF1, alpha_VBF2);
         qgl_VBF1 = AK4puppi_VBF[0].QGL();
         qgl_VBF2 = AK4puppi_VBF[1].QGL();
         QvsG_VBF1 = AK4puppi_VBF[0].QvsG();
         QvsG_VBF2 = AK4puppi_VBF[1].QvsG();
         pz_4jets = jetc1_p4.Pz() + jetc2_p4.Pz() + jetVBF1_p4.Pz() + jetVBF2_p4.Pz();
         TLorentzVector p4_4jets = jetc1_p4 + jetc2_p4 + jetVBF1_p4 + jetVBF2_p4;
         pt_norm = p4_4jets.Pt()/(jetc1_p4.Pt() + jetc2_p4.Pt() + jetVBF1_p4.Pt() + jetVBF2_p4.Pt());
         TLorentzVector Higgs_p4 = jetc1_p4 + jetc2_p4;
         DR_HiggsVBF1 = Higgs_p4.DeltaR(jetVBF1_p4);
         DR_HiggsVBF2 = Higgs_p4.DeltaR(jetVBF2_p4);
         Dphi_qq_cc = fabs(VBFsum.DeltaPhi(Higgs_p4));
         mCC= (jetc1_p4+jetc2_p4).M(); 
         tree->Fill();
       }
     }//end condition on HLT
     for(int i=0; i<NCUTS; i++){
       if(Eff_counter[i] == true) cutevt[i]++;
     }
     if(verbose){
      if(Eff_counter[4]==false && Eff_counter[5]==true){
          cout<<"AAAAAAAAAAAAAAAA"<<endl;
      }
     }
 } //end cycle on events

   if(Apply_jetVeto==true) vetoFile->Close();
   if(isMC!=0 && Apply_PUrew==true) PUrewFile->Close();
   fout->cd();
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
