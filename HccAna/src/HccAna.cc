// -*- C++ -*-
//
// Package:    HccAna
// Class:      HccAna
// 
///

// system include files
#include <memory>
#include <string>
#include <map>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <iomanip>
#include <set>

#define PI 3.14159

// user include files 
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TSpline.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "Math/VectorUtil.h"
#include "TClonesArray.h"
#include "TCanvas.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/limited/EDAnalyzerBase.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

//HTXS
#include "SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h"
//#include "SimDataFormats/HZZFiducial/interface/HZZFiducialVolume.h"

// PAT
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/Provenance/interface/Timestamp.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

//L1trigger
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/L1TGlobal/interface/GlobalExtBlk.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1TUtmTriggerMenu.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"

// Reco
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

//Helper
#include "Hcc/HccAna/interface/HccHelper.h"
//Muons
#include "Hcc/HccAna/interface/HccMuonAna.h"
#include "Hcc/HccAna/interface/HccMuonTree.h"
//Electrons
#include "Hcc/HccAna/interface/HccElectronTree.h"
//Photons
#include "Hcc/HccAna/interface/HccPhotonTree.h"
//Jets
#include "Hcc/HccAna/interface/HccJetTree.h"
//Final Leps
#include "Hcc/HccAna/interface/HccFinalLepTree.h"
//Sip
#include "Hcc/HccAna/interface/HccSipAna.h"
//PU
#include "Hcc/HccAna/interface/HccPileUp.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

//GEN
#include "Hcc/HccAna/interface/HccGENAna.h"
//VBF Jets
#include "Hcc/HccAna/interface/HccJets.h"

// Jet energy correction
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

// Jet flavour info
//#include "SimDataFormats/JetMatching/interface/JetFlavourInfo.h"

#include <vector>

// Kinematic Fit
#include "RecoParticleFlow/PFClusterTools/interface/PFEnergyResolution.h"

// EWK corrections
#include "Hcc/HccAna/interface/EwkCorrections.h"

// JEC related
//#include "PhysicsTools/PatAlgos/plugins/PATJetUpdater.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
 
//JER related
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

//BTag Calibration

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

//Muon MVA
//#include "MuonMVAReader/Reader/interface/MuonGBRForestReader.hpp"

// KalmanVertexFitter  
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
// Rochester Corrections
//#include "Hcc/KalmanMuonCalibrationsProducer/src/RoccoR.cc"

#include "RecoVertex/KalmanVertexFit/interface/SingleTrackVertexConstraint.h"

//
// class declaration
//
using namespace EwkCorrections;

class HccAna : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
    explicit HccAna(const edm::ParameterSet&);
    ~HccAna();
  
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    static bool sortByPt( const reco::GenParticle &p1, const reco::GenParticle &p2 ){ return (p1.pt() > p2.pt()); };
  
private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup& );
    virtual void endJob() ;
  
    virtual void beginRun(edm::Run const&, const edm::EventSetup& iSetup);
    virtual void endRun(edm::Run const&, edm::EventSetup const&);
    virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
    virtual void endLuminosityBlock(edm::LuminosityBlock const& lumiSeg,edm::EventSetup const& eSetup);
  
    //Helper Class
    HccHelper helper;
    //GEN
    HccGENAna genAna;
    //VBF
    HccJets jetHelper;
    //PU Reweighting
    edm::LumiReWeighting *lumiWeight;
    HccPileUp pileUp;
    //JES Uncertainties
    std::unique_ptr<JetCorrectionUncertainty> jecunc;
    // kfactors
    TSpline3 *kFactor_ggzz;
    std::vector<std::vector<float> > tableEwk;
    // data/MC scale factors
    TH2F *hElecScaleFac;
    TH2F *hElecScaleFac_Cracks;
    TH2F *hElecScaleFacGsf;
    TH2F *hElecScaleFacGsfLowET;
    TH2F *hMuScaleFac;
    TH2F *hMuScaleFacUnc;
    TH1D *h_pileup;
    TH1D *h_pileupUp;
    TH1D *h_pileupDn;
    std::vector<TH1F*> h_medians;
    TH2F *hbTagEffi;
    TH2F *hcTagEffi;
    TH2F *hudsgTagEffi;

    BTagCalibrationReader* reader;

    //Saved Events Trees
    TTree *passedEventsTree_All;

    void bookPassedEventTree(TString treeName, TTree *tree);
    /*void setTreeVariables( const edm::Event&, const edm::EventSetup&, 
                           std::vector<pat::Muon> selectedMuons, std::vector<pat::Electron> selectedElectrons, 
                           std::vector<pat::Muon> recoMuons, std::vector<pat::Electron> recoElectrons, 
                           std::vector<pat::Jet> goodJets, std::vector<float> goodJetQGTagger, 
                           std::vector<float> goodJetaxis2, std::vector<float> goodJetptD, std::vector<int> goodJetmult, 
                           std::vector<pat::Jet> selectedMergedJets,
                           std::map<unsigned int, TLorentzVector> selectedFsrMap);
    void setGENVariables(edm::Handle<reco::GenParticleCollection> prunedgenParticles,
                         edm::Handle<edm::View<pat::PackedGenParticle> > packedgenParticles,
                         edm::Handle<edm::View<reco::GenJet> > genJets);*/
		void setTreeVariables( const edm::Event&, const edm::EventSetup&,
                           //std::vector<pat::Jet> goodJets,// std::vector<float> goodJetQGTagger,
                           //std::vector<float> goodJetaxis2, std::vector<float> goodJetptD, std::vector<int> goodJetmult,
                           //std::vector<pat::Jet> selectedMergedJets,
                           edm::Handle<edm::View<pat::Jet> > AK4PuppiJets,
                           //edm::Handle<edm::View<pat::Jet> > AK4PuppiJetsJESUp,
                           //edm::Handle<edm::View<pat::Jet> > AK4PuppiJetsJESDown,
                           edm::Handle<edm::View<pat::Jet> > AK4PuppiJetsUncorr,
                           JetCorrectionUncertainty *uncAK4,
                           edm::Handle<edm::View<pat::Jet> > AK4PuppiJetsSmear,
                           edm::Handle<edm::View<pat::Jet> > AK4PuppiJetsSmearUp,
                           edm::Handle<edm::View<pat::Jet> > AK4PuppiJetsSmearDown,
                           edm::Handle<edm::View<pat::Jet> > AK8PuppiJets,
                           //const edm::TriggerNames trigNames,
                           //edm::Handle<pat::TriggerObjectStandAlone> triggerObjects,
                           //edm::Handle<std::vector<reco::PFJet>> hltjets,
                           //edm::Handle<edm::View<reco::PFJet>> hltjetsForBTag,
                           //edm::Handle<edm::View<reco::PFJet>> hltAK4PFJetsCorrected,
                           //edm::Handle<reco::JetTagCollection> pfJetTagCollectionPrticleNetprobc,
                           //edm::Handle<reco::JetTagCollection> pfJetTagCollectionPrticleNetprobb,
                           //edm::Handle<reco::JetTagCollection> pfJetTagCollectionPrticleNetprobuds,
                           //edm::Handle<reco::JetTagCollection> pfJetTagCollectionPrticleNetprobg,
                           //edm::Handle<reco::JetTagCollection> pfJetTagCollectionPrticleNetprobtauh,
                           //edm::Handle<BXVector<l1t::Jet> > bxvCaloJets,
                           //edm::Handle<BXVector<l1t::Muon> > bxvCaloMuons,
                           //edm::Handle<BXVector<l1t::EtSum> > bxvCaloHT,
                           std::vector<pat::Muon> AllMuons, std::vector<pat::Electron> AllElectrons, const reco::Vertex *ver);
    void setGENVariables(edm::Handle<reco::GenParticleCollection> prunedgenParticles,
                         edm::Handle<edm::View<pat::PackedGenParticle> > packedgenParticles,
                         edm::Handle<edm::View<reco::GenJet> > genJets);


    // -------------------------
    // RECO level information
    // -------------------------

    // Event Variables
    ULong64_t Run, Event, LumiSect, puN;
    int nVtx, nInt;
    int finalState;
    std::string triggersPassed;
    bool passedTrig, passedFullSelection, passedZ4lSelection, passedQCDcut, passedZqqSelection;

    std::vector<string>  Trigger_l1name;
    std::vector<int> Trigger_l1decision;
    std::vector<string>  Trigger_hltname;
    std::vector<int> Trigger_hltdecision;
    
    float PV_x, PV_y, PV_z; 
    float BS_x, BS_y, BS_z; 
    float BS_xErr, BS_yErr, BS_zErr; 
    float BeamWidth_x, BeamWidth_y;
    float BeamWidth_xErr, BeamWidth_yErr;


    // Event Weights
    float genWeight, pileupWeight, pileupWeightUp, pileupWeightDn, dataMCWeight, eventWeight, prefiringWeight;
    float k_qqZZ_qcd_dPhi, k_qqZZ_qcd_M, k_qqZZ_qcd_Pt, k_qqZZ_ewk;
    // pdf weights                                                                   
    vector<float> qcdWeights;
    vector<float> nnloWeights;
    vector<float> pdfWeights;
    int posNNPDF;
    float pdfRMSup, pdfRMSdown, pdfENVup, pdfENVdown;

    // lepton variables
    //vector<double> ALLlep_pt; vector<double> ALLlep_eta; vector<double> ALLlep_phi; vector<double> ALLlep_mass; vector<int> ALLlep_id;
    vector<double> Ele_pt; vector<double> Ele_eta; vector<double> Ele_phi; vector<double> Ele_mass; vector<double> Ele_dxy; vector<double> Ele_dz; vector<int> Ele_id; vector<double> Ele_hcalIso; vector<double> Ele_ecalIso; vector<double> Ele_trackIso; vector<bool> Ele_isEB; vector<double> Ele_IsoCal; vector<double> Ele_03_Neutral; 
/*vector<double> Ele_PF_Iso_R04;*/ vector<bool> Ele_isPassID;
    vector<double> Muon_pt; vector<double> Muon_eta; vector<double> Muon_phi; vector<double> Muon_mass; vector<double> Muon_dxy; vector<double> Muon_dz; vector<int> Muon_id; vector<bool> Muon_isPF; vector<double> Muon_PF_Iso_R04; vector<bool> Muon_PassLooseID;
    vector<double> AK4lep_pt; vector<double> AK4lep_eta; vector<double> AK4lep_phi; vector<double> AK4lep_mass; vector<int> AK4lep_id;

    // MET
    float met; float met_pt; float met_phi;
    //float met_jesup, met_phi_jesup, met_jesdn, met_phi_jesdn;
    //float met_uncenup, met_phi_uncenup, met_uncendn, met_phi_uncendn;

    //L1 HT
    float L1ht;

    //hlt jets for B Tag
    vector<double> hltjetForBTag_pt;
    vector<double> hltjetForBTag_eta;
    vector<double> hltjetForBTag_phi;
    vector<double> hltjetForBTag_mass;
    vector<float> hltParticleNetONNXJetTags_probb, hltParticleNetONNXJetTags_probc,hltParticleNetONNXJetTags_probuds, hltParticleNetONNXJetTags_probg, hltParticleNetONNXJetTags_probtauh;    
		
    // HLT jets hltAK4PFJetsCorrected
		
    vector<double> hltAK4PFJetsCorrected_pt;
    vector<double> hltAK4PFJetsCorrected_eta;
    vector<double> hltAK4PFJetsCorrected_phi;
    vector<double> hltAK4PFJetsCorrected_mass;
   
    //HLT jets for turn on curves
    //HLT_PFJet80
    vector<double> HLTJet80_pt, HLTJet80_eta, HLTJet80_phi;
    vector<double> HLTJet60_pt, HLTJet60_eta, HLTJet60_phi;
    vector<double> HLTJet80_MatchedCalo_pt, HLTJet80_MatchedCalo_eta, HLTJet80_MatchedCalo_phi;  //this hlt object does not have the cut on pt=80
    vector<double> HLTJet60_MatchedCalo_pt, HLTJet60_MatchedCalo_eta, HLTJet60_MatchedCalo_phi;  //this hlt object does not have the cut on pt=80
    vector<double> HLTAK4PFJetLoose_pt, HLTAK4PFJetLoose_eta, HLTAK4PFJetLoose_phi;
    vector<double> HLTAK4PFJetTight_pt, HLTAK4PFJetTight_eta, HLTAK4PFJetTight_phi;
    vector<double> HLTAK4PFJet_pt, HLTAK4PFJet_eta, HLTAK4PFJet_phi;
 
    // Puppi AK4jets with ParticleNet taggers

    vector<double> AK4PuppiJets_pt;
    vector<double> AK4PuppiJetsJESUp_pt;
    vector<double> AK4PuppiJetsJESDown_pt;
    vector<double> AK4PuppiJetsSmear_pt;
    vector<double> AK4PuppiJetsSmearUp_pt;
    vector<double> AK4PuppiJetsSmearDown_pt;
    vector<double> AK4PuppiJets_uncorrpt;
    vector<double> AK4PuppiJets_eta;
    vector<double> AK4PuppiJets_phi;
    vector<double> AK4PuppiJets_mass;
    vector<int> AK4PuppiJets_hadrFlav;
    vector<int> AK4PuppiJets_partFlav;

    vector<float> jet_pfParticleNetAK4JetTags_probb, jet_pfParticleNetAK4JetTags_probc, jet_pfParticleNetAK4JetTags_probuds,jet_pfParticleNetAK4JetTags_probg, jet_pfParticleNetAK4JetTags_probtauh;  
    vector<float> jet_pfParticleNetAK4JetTags_CvsB, jet_pfParticleNetAK4JetTags_CvsL, jet_pfParticleNetAK4JetTags_CvsAll,jet_pfParticleNetAK4JetTags_BvsC, jet_pfParticleNetAK4JetTags_BvsL, jet_pfParticleNetAK4JetTags_BvsAll, jet_pfParticleNetAK4JetTags_QvsG; 
    vector<float> jet_PNetRegPtRawCorr, jet_PNetRegPtRawCorrNeutrino, jet_PNetRegPtRawRes; 


    vector<float> jet_pfDeepJetAK4JetTags_probb, jet_pfDeepJetAK4JetTags_probbb, jet_pfDeepJetAK4JetTags_problepb, jet_pfDeepJetAK4JetTags_probc, jet_pfDeepJetAK4JetTags_probuds,jet_pfDeepJetAK4JetTags_probg; 

    vector<float> jet_pfDeepCSVAK4JetTags_probb, jet_pfDeepCSVAK4JetTags_probbb, jet_pfDeepCSVAK4JetTags_probc, jet_pfDeepCSVAK4JetTags_probudsg; 
    // Puppi AK8jets with ParticleNet(-MD) and DeepDoubleX taggers

    int leadingAK8_pt_idx;
    int subleadingAK8_pt_idx;
	
    vector<double> AK8PuppiJets_pt;
    vector<double> AK8PuppiJets_eta;
	vector<double> AK8PuppiJets_phi;
	vector<double> AK8PuppiJets_mass;
    vector<double> AK8PuppiJets_softdropmass;
	
	vector<float> jet_pfParticleNetJetTags_probZbb, jet_pfParticleNetJetTags_probZcc, jet_pfParticleNetJetTags_probZqq, jet_pfParticleNetJetTags_probQCDbb, jet_pfParticleNetJetTags_probQCDcc, jet_pfParticleNetJetTags_probQCDb, jet_pfParticleNetJetTags_probQCDc, jet_pfParticleNetJetTags_probQCDothers, jet_pfParticleNetJetTags_probHbb, jet_pfParticleNetJetTags_probHcc, jet_pfParticleNetJetTags_probHqqqq;  
		
	vector<float> jet_pfMassDecorrelatedParticleNetJetTags_probXbb, jet_pfMassDecorrelatedParticleNetJetTags_probXcc, jet_pfMassDecorrelatedParticleNetJetTags_probXqq, jet_pfMassDecorrelatedParticleNetJetTags_probQCDbb, jet_pfMassDecorrelatedParticleNetJetTags_probQCDcc, jet_pfMassDecorrelatedParticleNetJetTags_probQCDb, jet_pfMassDecorrelatedParticleNetJetTags_probQCDc, jet_pfMassDecorrelatedParticleNetJetTags_probQCDothers;
    vector<float> jet_pfMassIndependentDeepDoubleBvLV2JetTags_probHbb, jet_pfMassIndependentDeepDoubleCvLV2JetTags_probHcc, jet_pfMassIndependentDeepDoubleCvBV2JetTags_probHcc;

    // Jets
    vector<int>    jet_iscleanH4l;
    int jet1index, jet2index;

	vector<double> L1jet_pt; vector<double> L1jet_eta; vector<double> L1jet_phi; vector<double> L1jet_mass;
    vector<double> L1muon_pt; vector<double> L1muon_eta; vector<double> L1muon_phi; vector<double> L1muon_mass;
	vector<int> L1muon_qual;	

    // Event Category
    int EventCat;

    // -------------------------
    // GEN level information
    // -------------------------

    //Event variables
    int GENfinalState;

    // Jets
    vector<double> GENjet_pt; vector<double> GENjet_eta; vector<double> GENjet_phi; vector<double> GENjet_mass; 
    vector<double> quark_pt; vector<double> quark_eta; vector<double> quark_phi; vector<int> quark_flavour; vector<bool> quark_VBF;
    float Z_pt; float Z_eta; float Z_phi; float Z_mass;
    int GENnjets_pt30_eta4p7; float GENpt_leadingjet_pt30_eta4p7; 
    int GENnjets_pt30_eta2p5; float GENpt_leadingjet_pt30_eta2p5; 
    float GENabsrapidity_leadingjet_pt30_eta4p7; float GENabsdeltarapidity_hleadingjet_pt30_eta4p7;
    int lheNb, lheNj, nGenStatus2bHad;

    int n_jets=0;
    vector<float> hltjetForBTag_pt_float, hltjetForBTag_eta_float, hltjetForBTag_phi_float, hltjetForBTag_mass_float;
    vector<float> hltAK4PFJetsCorrected_pt_float, hltAK4PFJetsCorrected_eta_float, hltAK4PFJetsCorrected_phi_float, hltAK4PFJetsCorrected_mass_float;
    vector<float> hltAK8PFJetsCorrected_pt_float, hltAK8PFJetsCorrected_eta_float, hltAK8PFJetsCorrected_phi_float, hltAK8PFJetsCorrected_mass_float;

    int n_GENjets=0;
    vector<float> GENjet_pt_float, GENjet_eta_float;
    vector<float> GENjet_phi_float, GENjet_mass_float;
	vector<float> quark_pt_float, quark_eta_float, quark_phi_float;
    vector<float> L1jet_pt_float, L1jet_eta_float, L1jet_phi_float, L1jet_mass_float;
    vector<float> L1muon_pt_float, L1muon_eta_float, L1muon_phi_float, L1muon_mass_float;

    vector<float> AK4PuppiJets_pt_float;
    vector<float> AK4PuppiJetsJESUp_pt_float;
    vector<float> AK4PuppiJetsJESDown_pt_float;
    vector<float> AK4PuppiJetsSmear_pt_float;
    vector<float> AK4PuppiJetsSmearUp_pt_float;
    vector<float> AK4PuppiJetsSmearDown_pt_float;
    vector<float> AK4PuppiJets_uncorrpt_float;
    vector<float> AK4PuppiJets_eta_float;
    vector<float> AK4PuppiJets_phi_float;
    vector<float> AK4PuppiJets_mass_float;
    vector<float> AK4PuppiJets_qgl;
	
    vector<float> HLTJet80_pt_float;
    vector<float> HLTJet80_eta_float;
    vector<float> HLTJet80_phi_float;

    vector<float> HLTJet60_pt_float;
    vector<float> HLTJet60_eta_float;
    vector<float> HLTJet60_phi_float;

    vector<float> HLTJet80_MatchedCalo_pt_float;
    vector<float> HLTJet80_MatchedCalo_eta_float;
    vector<float> HLTJet80_MatchedCalo_phi_float;

    vector<float> HLTJet60_MatchedCalo_pt_float;
    vector<float> HLTJet60_MatchedCalo_eta_float;
    vector<float> HLTJet60_MatchedCalo_phi_float;

    vector<float> HLTAK4PFJetLoose_pt_float;
    vector<float> HLTAK4PFJetLoose_eta_float;
    vector<float> HLTAK4PFJetLoose_phi_float;

    vector<float> HLTAK4PFJetTight_pt_float;
    vector<float> HLTAK4PFJetTight_eta_float;
    vector<float> HLTAK4PFJetTight_phi_float;

    vector<float> HLTAK4PFJet_pt_float;
    vector<float> HLTAK4PFJet_eta_float;
    vector<float> HLTAK4PFJet_phi_float;

	vector<float> AK8PuppiJets_pt_float;
    vector<float> AK8PuppiJets_eta_float;
    vector<float> AK8PuppiJets_phi_float;
    vector<float> AK8PuppiJets_mass_float;

    // Global Variables but not stored in the tree
    //vector<double> lep_ptreco;
    //vector<int> lep_ptid; vector<int> lep_ptindex;
    vector<pat::Muon> recoMuons; vector<pat::Electron> recoElectrons; vector<pat::Electron> recoElectronsUnS; 
    /*vector<pat::Tau> recoTaus; vector<pat::Photon> recoPhotons;
    vector<pat::PFParticle> fsrPhotons; 
    TLorentzVector HVec, HVecNoFSR, Z1Vec, Z2Vec;
    TLorentzVector GENZ1Vec, GENZ2Vec;
    bool foundHiggsCandidate; bool firstEntry;*/
    float jet1pt, jet2pt;
		bool firstEntry;

    // hist container
    std::map<std::string,TH1F*> histContainer_;

    //Input edm
    edm::EDGetTokenT<edm::View<pat::Electron> > elecSrc_;
    edm::EDGetTokenT<edm::View<pat::Electron> > elecUnSSrc_;
    edm::EDGetTokenT<edm::View<pat::Muon> > muonSrc_;
    //edm::EDGetTokenT<edm::View<pat::Photon> > photonSrc_;
    //edm::EDGetTokenT<edm::View<pat::Jet> > jetSrc_;
    edm::EDGetTokenT<edm::View<pat::Jet> > AK4PuppiJetSrc_;
    //edm::EDGetTokenT<edm::View<pat::Jet> > AK4PuppiJetJESUpSrc_;
    //edm::EDGetTokenT<edm::View<pat::Jet> > AK4PuppiJetJESDownSrc_;
    edm::EDGetTokenT<edm::View<pat::Jet> > AK4PuppiJetUncorrSrc_;
    edm::EDGetTokenT<edm::View<pat::Jet> > AK4PuppiJetSmearSrc_;
    edm::EDGetTokenT<edm::View<pat::Jet> > AK4PuppiJetSmearUpSrc_;
    edm::EDGetTokenT<edm::View<pat::Jet> > AK4PuppiJetSmearDownSrc_;
    edm::EDGetTokenT<edm::View<pat::Jet> > AK8PuppiJetSrc_;
    //edm::EDGetTokenT<BXVector<l1t::Jet>> bxvCaloJetSrc_;
    //edm::EDGetTokenT<edm::View<reco::PFJet>> hltPFJetForBtagSrc_;
    //edm::EDGetTokenT<edm::View<reco::PFJet>> hltAK4PFJetsCorrectedSrc_;
    //edm::EDGetTokenT<reco::JetTagCollection> pfJetTagCollectionParticleNetprobcSrc_;  //value map for Particle Net tagger at hlt
    //edm::EDGetTokenT<reco::JetTagCollection> pfJetTagCollectionParticleNetprobbSrc_;  //value map for Particle Net tagger at hlt
    //edm::EDGetTokenT<reco::JetTagCollection> pfJetTagCollectionParticleNetprobudsSrc_;  //value map for Particle Net tagger at hlt
    //edm::EDGetTokenT<reco::JetTagCollection> pfJetTagCollectionParticleNetprobgSrc_;  //value map for Particle Net tagger at hlt
    //edm::EDGetTokenT<reco::JetTagCollection> pfJetTagCollectionParticleNetprobtauhSrc_;  //value map for Particle Net tagger at hlt
    //edm::EDGetTokenT<BXVector<l1t::Muon>> bxvCaloMuonSrc_;
    //edm::EDGetTokenT<BXVector<l1t::EtSum>> bxvCaloHTSrc_;
    edm::EDGetTokenT<edm::ValueMap<float> > qgTagSrc_;
    edm::EDGetTokenT<edm::ValueMap<float> > axis2Src_;
    edm::EDGetTokenT<edm::ValueMap<int> > multSrc_;
    edm::EDGetTokenT<edm::ValueMap<float> > ptDSrc_;
    //edm::EDGetTokenT<edm::View<pat::Jet> > mergedjetSrc_;
    edm::EDGetTokenT<edm::View<pat::MET> > metSrc_;
    //edm::InputTag triggerSrc_;
    edm::EDGetTokenT<edm::TriggerResults> triggerSrc_;
    //edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
    edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObjects_;
    edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;
    edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;
    edm::EDGetTokenT<std::vector<reco::Conversion> > conversionSrc_;
    edm::EDGetTokenT<double> muRhoSrc_;
    edm::EDGetTokenT<double> elRhoSrc_;
    edm::EDGetTokenT<double> rhoSrcSUS_;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupSrc_;
    edm::EDGetTokenT<pat::PackedCandidateCollection> pfCandsSrc_;
    edm::EDGetTokenT<edm::View<pat::PFParticle> > fsrPhotonsSrc_;
    edm::EDGetTokenT<reco::GenParticleCollection> prunedgenParticlesSrc_;
    edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedgenParticlesSrc_;
    edm::EDGetTokenT<edm::View<reco::GenJet> > genJetsSrc_;
    edm::EDGetTokenT<GenEventInfoProduct> generatorSrc_;
    edm::EDGetTokenT<LHEEventProduct> lheInfoSrc_;
    edm::EDGetTokenT<LHERunInfoProduct> lheRunInfoToken_;
    edm::EDGetTokenT<HTXS::HiggsClassification> htxsSrc_;
    //edm::EDGetTokenT<HZZFid::FiducialSummary> fidRivetSrc_;
    edm::EDGetTokenT< double > prefweight_token_;
    //L1 Trigger
    //edm::EDGetToken algTok_;
    //edm::EDGetTokenT<GlobalAlgBlkBxCollection> algInputTag_;
    //l1t::L1TGlobalUtil* gtUtil_;


    // Configuration
    const float Zmass;
    float mZ1Low, mZ2Low, mZ1High, mZ2High, m4lLowCut;
    float jetpt_cut, jeteta_cut;
    std::string elecID;
    bool isMC, isSignal;
    bool isHcc, isZqq, isZcc, isZbb;
    float mH;
    float crossSection;
    bool weightEvents;
    float isoCutEl, isoCutMu; 
    double isoConeSizeEl, isoConeSizeMu;
    float sip3dCut, leadingPtCut, subleadingPtCut;
    float genIsoCutEl, genIsoCutMu;
    double genIsoConeSizeEl, genIsoConeSizeMu;
    float _elecPtCut, _muPtCut, _tauPtCut, _phoPtCut;
    float BTagCut;
    bool reweightForPU;
    std::string PUVersion;
    bool doFsrRecovery, GENbestM4l;
    bool doPUJetID;
    int jetIDLevel;
    bool doJER;
    bool doJEC;
    bool doRefit;
    bool doTriggerMatching;
    bool checkOnlySingle;
    std::vector<std::string> triggerList;
    int skimLooseLeptons, skimTightLeptons;
    bool verbose;

    int year;///use to choose Muon BDT
    bool isCode4l;

   std::string JECUncFileAK4;
   std::vector<JetCorrectionUncertainty*> vsrc ;

edm::ESGetToken<JetCorrectorParametersCollection, JetCorrectionsRecord> mPayloadToken;

std::string res_pt_config;
std::string res_phi_config;
std::string res_sf_config;

    // register to the TFileService
    edm::Service<TFileService> fs;

    // Counters
    float nEventsTotal;
    float sumWeightsTotal;
    float sumWeightsTotalPU;

    // JER
    JME::JetResolution resolution_pt, resolution_phi;
    JME::JetResolutionScaleFactor resolution_sf;

    string EleBDT_name_161718;
    string heepID_name_161718;

    //JEC uncertainty files
    JetCorrectorParameters *pAK4 = new JetCorrectorParameters(JECUncFileAK4.c_str(), "Total") ; 
    JetCorrectionUncertainty *uncAK4 = new JetCorrectionUncertainty(*pAK4);
};


HccAna::HccAna(const edm::ParameterSet& iConfig) :
    histContainer_(),
    elecSrc_(consumes<edm::View<pat::Electron> >(iConfig.getUntrackedParameter<edm::InputTag>("electronSrc"))),
    //elecSrc_(consumes<edm::View<pat::Electron> >(iConfig.getUntrackedParameter<edm::InputTag>("electronUnSSrc"))),
    elecUnSSrc_(consumes<edm::View<pat::Electron> >(iConfig.getUntrackedParameter<edm::InputTag>("electronUnSSrc"))),
    muonSrc_(consumes<edm::View<pat::Muon> >(iConfig.getUntrackedParameter<edm::InputTag>("muonSrc"))),
    //photonSrc_(consumes<edm::View<pat::Photon> >(iConfig.getUntrackedParameter<edm::InputTag>("photonSrc"))),
    //jetSrc_(consumes<edm::View<pat::Jet> >(iConfig.getUntrackedParameter<edm::InputTag>("jetSrc"))),
    AK4PuppiJetSrc_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("AK4PuppiJetSrc"))),
    //AK4PuppiJetJESUpSrc_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("AK4PuppiJetJESUpSrc"))),
    //AK4PuppiJetJESDownSrc_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("AK4PuppiJetJESDownSrc"))),
    AK4PuppiJetUncorrSrc_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("AK4PuppiJetUncorrSrc"))),
    AK4PuppiJetSmearSrc_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("AK4PuppiJetSmearSrc"))),
    AK4PuppiJetSmearUpSrc_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("AK4PuppiJetSmearUpSrc"))),
    AK4PuppiJetSmearDownSrc_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("AK4PuppiJetSmearDownSrc"))),
    //AK4PuppiJetSrc_(consumes<edm::View<pat::Jet> >(iConfig.getUntrackedParameter<edm::InputTag>("AK4PuppiJetSrc"))),
	AK8PuppiJetSrc_(consumes<edm::View<pat::Jet> >(iConfig.getUntrackedParameter<edm::InputTag>("AK8PuppiJetSrc"))),
    //bxvCaloJetSrc_(consumes<BXVector<l1t::Jet>>(iConfig.getParameter<edm::InputTag>("bxvCaloJetSrc"))),
    //hltPFJetForBtagSrc_(consumes<edm::View<reco::PFJet>>(iConfig.getParameter<edm::InputTag>("hltPFJetForBtagSrc"))),
    //hltAK4PFJetsCorrectedSrc_(consumes<edm::View<reco::PFJet>>(iConfig.getParameter<edm::InputTag>("hltAK4PFJetsCorrectedSrc"))),
    //pfJetTagCollectionParticleNetprobcSrc_(consumes(iConfig.getParameter<edm::InputTag>("pfJetTagCollectionParticleNetprobcSrc"))),
    //pfJetTagCollectionParticleNetprobbSrc_(consumes(iConfig.getParameter<edm::InputTag>("pfJetTagCollectionParticleNetprobbSrc"))),
    //pfJetTagCollectionParticleNetprobudsSrc_(consumes(iConfig.getParameter<edm::InputTag>("pfJetTagCollectionParticleNetprobudsSrc"))),
    //pfJetTagCollectionParticleNetprobgSrc_(consumes(iConfig.getParameter<edm::InputTag>("pfJetTagCollectionParticleNetprobgSrc"))),
    //pfJetTagCollectionParticleNetprobtauhSrc_(consumes(iConfig.getParameter<edm::InputTag>("pfJetTagCollectionParticleNetprobtauhSrc"))),
    //bxvCaloMuonSrc_(consumes<BXVector<l1t::Muon>>(iConfig.getParameter<edm::InputTag>("bxvCaloMuonSrc"))),
    //bxvCaloHTSrc_(consumes<BXVector<l1t::EtSum>>(iConfig.getParameter<edm::InputTag>("bxvCaloHTSrc"))),
    //qgTagSrc_(consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "qgLikelihood"))),
    //axis2Src_(consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "axis2"))),
    //multSrc_(consumes<edm::ValueMap<int>>(edm::InputTag("QGTagger", "mult"))),
    //ptDSrc_(consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "ptD"))),
    //mergedjetSrc_(consumes<edm::View<pat::Jet> >(iConfig.getUntrackedParameter<edm::InputTag>("mergedjetSrc"))),
    metSrc_(consumes<edm::View<pat::MET> >(iConfig.getUntrackedParameter<edm::InputTag>("metSrc"))),
    triggerSrc_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerSrc"))),
    //triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone> >(iConfig.getParameter<edm::InputTag>("objects")));
    triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
    //triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
    vertexSrc_(consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("vertexSrc"))),
    beamSpotSrc_(consumes<reco::BeamSpot>(iConfig.getUntrackedParameter<edm::InputTag>("beamSpotSrc"))),
    conversionSrc_(consumes<std::vector<reco::Conversion> >(iConfig.getUntrackedParameter<edm::InputTag>("conversionSrc"))),
    muRhoSrc_(consumes<double>(iConfig.getUntrackedParameter<edm::InputTag>("muRhoSrc"))),
    elRhoSrc_(consumes<double>(iConfig.getUntrackedParameter<edm::InputTag>("elRhoSrc"))),
    rhoSrcSUS_(consumes<double>(iConfig.getUntrackedParameter<edm::InputTag>("rhoSrcSUS"))),
    pileupSrc_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getUntrackedParameter<edm::InputTag>("pileupSrc"))),
    pfCandsSrc_(consumes<pat::PackedCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("pfCandsSrc"))),
    fsrPhotonsSrc_(consumes<edm::View<pat::PFParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("fsrPhotonsSrc"))),
    prunedgenParticlesSrc_(consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("prunedgenParticlesSrc"))),
    packedgenParticlesSrc_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("packedgenParticlesSrc"))),
    genJetsSrc_(consumes<edm::View<reco::GenJet> >(iConfig.getUntrackedParameter<edm::InputTag>("genJetsSrc"))),
    generatorSrc_(consumes<GenEventInfoProduct>(iConfig.getUntrackedParameter<edm::InputTag>("generatorSrc"))),
    lheInfoSrc_(consumes<LHEEventProduct>(iConfig.getUntrackedParameter<edm::InputTag>("lheInfoSrc"))),
    lheRunInfoToken_(consumes<LHERunInfoProduct,edm::InRun>(edm::InputTag("externalLHEProducer",""))),
    htxsSrc_(consumes<HTXS::HiggsClassification>(edm::InputTag("rivetProducerHTXS","HiggsClassification"))),
    //prefweight_token_(consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProb"))),
    //fidRivetSrc_(consumes<HZZFid::FiducialSummary>(edm::InputTag("rivetProducerHZZFid","FiducialSummary"))),
    //L1 Trigger
    //algTok_(consumes<BXVector<GlobalAlgBlk>>(iConfig.getParameter<edm::InputTag>("algInputTag"))),
    //algInputTag_(consumes<GlobalAlgBlkBxCollection>(iConfig.getParameter<edm::InputTag>("algInputTag"))),
    //gtUtil_(new l1t::L1TGlobalUtil(iConfig, consumesCollector(), *this, iConfig.getParameter<edm::InputTag>("algInputTag"), iConfig.getParameter<edm::InputTag>("algInputTag"), l1t::UseEventSetupIn::RunAndEvent)),
    Zmass(91.1876),
    mZ1Low(iConfig.getUntrackedParameter<double>("mZ1Low",40.0)),
    mZ2Low(iConfig.getUntrackedParameter<double>("mZ2Low",12.0)), // was 12
    mZ1High(iConfig.getUntrackedParameter<double>("mZ1High",120.0)),
    mZ2High(iConfig.getUntrackedParameter<double>("mZ2High",120.0)),
    m4lLowCut(iConfig.getUntrackedParameter<double>("m4lLowCut",70.0)),
//     m4lLowCut(iConfig.getUntrackedParameter<double>("m4lLowCut",0.0)),
    jetpt_cut(iConfig.getUntrackedParameter<double>("jetpt_cut",10.0)),
    jeteta_cut(iConfig.getUntrackedParameter<double>("eta_cut",47)),
    elecID(iConfig.getUntrackedParameter<std::string>("elecID","NonTrig")),
    isMC(iConfig.getUntrackedParameter<bool>("isMC",true)),
    isSignal(iConfig.getUntrackedParameter<bool>("isSignal",false)),    
    isHcc(iConfig.getUntrackedParameter<bool>("isHcc")),
    isZqq(iConfig.getUntrackedParameter<bool>("isZqq")),
    isZcc(iConfig.getUntrackedParameter<bool>("isZcc")),
    isZbb(iConfig.getUntrackedParameter<bool>("isZbb")),
    mH(iConfig.getUntrackedParameter<double>("mH",0.0)),
    crossSection(iConfig.getUntrackedParameter<double>("CrossSection",1.0)),
    weightEvents(iConfig.getUntrackedParameter<bool>("weightEvents",false)),
    isoCutEl(iConfig.getUntrackedParameter<double>("isoCutEl",9999.0)),
    isoCutMu(iConfig.getUntrackedParameter<double>("isoCutMu",0.35)),/////ios is applied to new Muon BDT //previous 0.35///Qianying
    isoConeSizeEl(iConfig.getUntrackedParameter<double>("isoConeSizeEl",0.3)),
    isoConeSizeMu(iConfig.getUntrackedParameter<double>("isoConeSizeMu",0.3)),
    sip3dCut(iConfig.getUntrackedParameter<double>("sip3dCut",4)),
    leadingPtCut(iConfig.getUntrackedParameter<double>("leadingPtCut",20.0)),
    subleadingPtCut(iConfig.getUntrackedParameter<double>("subleadingPtCut",10.0)),
    genIsoCutEl(iConfig.getUntrackedParameter<double>("genIsoCutEl",0.35)), 
    genIsoCutMu(iConfig.getUntrackedParameter<double>("genIsoCutMu",0.35)), 
    genIsoConeSizeEl(iConfig.getUntrackedParameter<double>("genIsoConeSizeEl",0.3)), 
    genIsoConeSizeMu(iConfig.getUntrackedParameter<double>("genIsoConeSizeMu",0.3)), 
    _elecPtCut(iConfig.getUntrackedParameter<double>("_elecPtCut",0.0)),
    _muPtCut(iConfig.getUntrackedParameter<double>("_muPtCut",0.0)),
    _tauPtCut(iConfig.getUntrackedParameter<double>("_tauPtCut",20.0)),
    _phoPtCut(iConfig.getUntrackedParameter<double>("_phoPtCut",10.0)),
    BTagCut(iConfig.getUntrackedParameter<double>("BTagCut",0.4184)),/////2016: 0.6321; 2017: 0.4941; 2018: 0.4184
    reweightForPU(iConfig.getUntrackedParameter<bool>("reweightForPU",true)),
    PUVersion(iConfig.getUntrackedParameter<std::string>("PUVersion","Summer16_80X")),
    doFsrRecovery(iConfig.getUntrackedParameter<bool>("doFsrRecovery",true)),
    GENbestM4l(iConfig.getUntrackedParameter<bool>("GENbestM4l",false)),
    doPUJetID(iConfig.getUntrackedParameter<bool>("doPUJetID",true)),
    jetIDLevel(iConfig.getUntrackedParameter<int>("jetIDLevel",2)),
    doJER(iConfig.getUntrackedParameter<bool>("doJER",true)),
    doJEC(iConfig.getUntrackedParameter<bool>("doJEC",true)),
    doRefit(iConfig.getUntrackedParameter<bool>("doRefit",true)),
    doTriggerMatching(iConfig.getUntrackedParameter<bool>("doTriggerMatching",!isMC)),
    checkOnlySingle(iConfig.getUntrackedParameter<bool>("checkOnlySingle",false)),
    triggerList(iConfig.getUntrackedParameter<std::vector<std::string>>("triggerList")),
    skimLooseLeptons(iConfig.getUntrackedParameter<int>("skimLooseLeptons",2)),    
    skimTightLeptons(iConfig.getUntrackedParameter<int>("skimTightLeptons",2)),    
    verbose(iConfig.getUntrackedParameter<bool>("verbose",true)),
    year(iConfig.getUntrackedParameter<int>("year",2018)),
    ////for year put 
    // 20160 for pre VFP
    // 20165 for post VFP
    // 2017
    // 2018
    // to select correct training
    isCode4l(iConfig.getUntrackedParameter<bool>("isCode4l",true)),
    JECUncFileAK4(iConfig.getParameter<std::string>("JECUncFileAK4Src")),
    mPayloadToken    {esConsumes(edm::ESInputTag("", iConfig.getParameter<std::string>("payload")))}
{
  
    if(!isMC){reweightForPU = false;}

    nEventsTotal=0.0;
    sumWeightsTotal=0.0;
    sumWeightsTotalPU=0.0;
    histContainer_["NEVENTS"]=fs->make<TH1F>("nEvents","nEvents in Sample",2,0,2);
    histContainer_["SUMWEIGHTS"]=fs->make<TH1F>("sumWeights","sum Weights of Sample",2,0,2);
    histContainer_["SUMWEIGHTSPU"]=fs->make<TH1F>("sumWeightsPU","sum Weights and PU of Sample",2,0,2);
    histContainer_["NVTX"]=fs->make<TH1F>("nVtx","Number of Vertices",76,-0.5,75.5);
    histContainer_["NVTX_RW"]=fs->make<TH1F>("nVtx_ReWeighted","Number of Vertices",76,-0.5,75.5);
    histContainer_["NINTERACT"]=fs->make<TH1F>("nInteractions","Number of True Interactions",81,-0.5,80.5);
    histContainer_["NINTERACT_RW"]=fs->make<TH1F>("nInteraction_ReWeighted","Number of True Interactions",81,-0.5,80.5);

    passedEventsTree_All = new TTree("passedEvents","passedEvents");

    tableEwk = readFile_and_loadEwkTable("ZZBG");  //LUIGI 

	int YEAR = year - 2016 + 1;
	if(year == 20165) YEAR = 1;
	if(year == 20160) YEAR = 0;
    if(year == 2022 || year == 2023) YEAR = 3; //DA MODIFICARE
   
    //string elec_scalefac_Cracks_name_161718[3] = {"egammaEffi.txt_EGM2D_cracks.root", "egammaEffi.txt_EGM2D_Moriond2018v1_gap.root", "egammaEffi.txt_EGM2D_Moriond2019_v1_gap.root"};
    /*string elec_scalefac_Cracks_name_161718[4] = {"ElectronSF_UL2016preVFP_gap.root", "ElectronSF_UL2016postVFP_gap.root", "ElectronSF_UL2017_gap.root", "ElectronSF_UL2018_gap.root"};
    edm::FileInPath elec_scalefacFileInPathCracks(("Hcc/HccAna/data/"+elec_scalefac_Cracks_name_161718[YEAR]).c_str());
    TFile *fElecScalFacCracks = TFile::Open(elec_scalefacFileInPathCracks.fullPath().c_str());
    hElecScaleFac_Cracks = (TH2F*)fElecScalFacCracks->Get("EGamma_SF2D");    
    //string elec_scalefac_name_161718[3] = {"egammaEffi.txt_EGM2D.root", "egammaEffi.txt_EGM2D_Moriond2018v1.root", "egammaEffi.txt_EGM2D_Moriond2019_v1.root"};
    string elec_scalefac_name_161718[4] = {"ElectronSF_UL2016preVFP_nogap.root", "ElectronSF_UL2016postVFP_nogap.root", "ElectronSF_UL2017_nogap.root", "ElectronSF_UL2018_nogap.root"};
    edm::FileInPath elec_scalefacFileInPath(("Hcc/HccAna/data/"+elec_scalefac_name_161718[YEAR]).c_str());
    TFile *fElecScalFac = TFile::Open(elec_scalefacFileInPath.fullPath().c_str());
    hElecScaleFac = (TH2F*)fElecScalFac->Get("EGamma_SF2D");    

    //string elec_Gsfscalefac_name_161718[3] = {"egammaEffi.txt_EGM2D_GSF.root", "egammaEffi.txt_EGM2D_Moriond2018v1_runBCDEF_passingRECO.root", "Ele_Reco_2018.root"};//was previous;
    string elec_Gsfscalefac_name_161718[4] = {"egammaEffi_ptAbove20.txt_EGM2D_UL2016preVFP.root", "egammaEffi_ptAbove20.txt_EGM2D_UL2016postVFP.root", "egammaEffi_ptAbove20.txt_EGM2D_UL2017.root", "egammaEffi_ptAbove20.txt_EGM2D_UL2018.root"};
    edm::FileInPath elec_GsfscalefacFileInPath(("Hcc/HccAna/data/"+elec_Gsfscalefac_name_161718[YEAR]).c_str());
    TFile *fElecScalFacGsf = TFile::Open(elec_GsfscalefacFileInPath.fullPath().c_str());
    hElecScaleFacGsf = (TH2F*)fElecScalFacGsf->Get("EGamma_SF2D");

    //string elec_GsfLowETscalefac_name_161718[3]= {"", "egammaEffi.txt_EGM2D_Moriond2018v1_runBCDEF_passingRECO_lowEt.root", "Ele_Reco_LowEt_2018.root"};//was previous
    string elec_GsfLowETscalefac_name_161718[4]= {"egammaEffi_ptBelow20.txt_EGM2D_UL2016preVFP.root", "egammaEffi_ptBelow20.txt_EGM2D_UL2016postVFP.root", "egammaEffi_ptBelow20.txt_EGM2D_UL2017.root", "egammaEffi_ptBelow20.txt_EGM2D_UL2018.root"};
    edm::FileInPath elec_GsfLowETscalefacFileInPath(("Hcc/HccAna/data/"+elec_GsfLowETscalefac_name_161718[YEAR]).c_str());
    TFile *fElecScalFacGsfLowET = TFile::Open(elec_GsfLowETscalefacFileInPath.fullPath().c_str());
    hElecScaleFacGsfLowET = (TH2F*)fElecScalFacGsfLowET->Get("EGamma_SF2D");

    //string mu_scalefac_name_161718[3] = {"final_HZZ_Moriond17Preliminary_v4.root", "ScaleFactors_mu_Moriond2018_final.root", "final_HZZ_muon_SF_2018RunA2D_ER_2702.root"};//was previous; 
//         string mu_scalefac_name_161718[3] = {"final_HZZ_SF_2016_legacy_mupogsysts.root", "final_HZZ_SF_2017_rereco_mupogsysts_3010.root", "final_HZZ_SF_2018_rereco_mupogsysts_3010.root"};
//         string mu_scalefac_name_161718[4] = {"final_HZZ_muon_SF_2016RunB2H_legacy_newLoose_newIso_paper.root", "final_HZZ_muon_SF_2016RunB2H_legacy_newLoose_newIso_paper.root", "final_HZZ_muon_SF_2017_newLooseIso_mupogSysts_paper.root", "final_HZZ_muon_SF_2018RunA2D_ER_newLoose_newIso_paper.root"};
        string mu_scalefac_name_161718[4] = {"final_HZZ_SF_2016UL_mupogsysts_newLoose.root","final_HZZ_SF_2016UL_mupogsysts_newLoose.root","final_HZZ_SF_2017UL_mupogsysts_newLoose.root","final_HZZ_SF_2018UL_mupogsysts_newLoose.root"};
    edm::FileInPath mu_scalefacFileInPath(("Hcc/HccAna/data/"+mu_scalefac_name_161718[YEAR]).c_str());
    TFile *fMuScalFac = TFile::Open(mu_scalefacFileInPath.fullPath().c_str());
    hMuScaleFac = (TH2F*)fMuScalFac->Get("FINAL");
    hMuScaleFacUnc = (TH2F*)fMuScalFac->Get("ERROR");*/

    //string pileup_name_161718[3] = {"puWeightsMoriond17_v2.root", "puWeightsMoriond18.root", "pu_weights_2018.root"};///was previous
//    string pileup_name_161718[3] = {"pu_weights_2016.root", "pu_weights_2017.root", "pu_weights_2018.root"};
    string pileup_name_161718[4] = {"pileup_UL_2016_1plusShift.root", "pileup_UL_2016_1plusShift.root", "pileup_UL_2017_1plusShift.root", "pileup_UL_2018_1plusShift.root"};
    edm::FileInPath pileup_FileInPath(("Hcc/HccAna/data/"+pileup_name_161718[YEAR]).c_str());
    TFile *f_pileup = TFile::Open(pileup_FileInPath.fullPath().c_str());
    h_pileup = (TH1D*)f_pileup->Get("weights");
    h_pileupUp = (TH1D*)f_pileup->Get("weights_varUp");
    h_pileupDn = (TH1D*)f_pileup->Get("weights_varDn");



    /*string bTagEffi_name_161718[4] = {"bTagEfficiencies_2016.root", "bTagEfficiencies_2016.root", "bTagEfficiencies_2017.root", "bTagEfficiencies_2018.root"};
    edm::FileInPath BTagEffiInPath(("Hcc/HccAna/data/"+bTagEffi_name_161718[YEAR]).c_str());
    TFile *fbTagEffi = TFile::Open(BTagEffiInPath.fullPath().c_str());
    hbTagEffi = (TH2F*)fbTagEffi->Get("eff_b_M_ALL");
    hcTagEffi = (TH2F*)fbTagEffi->Get("eff_c_M_ALL");
    hudsgTagEffi = (TH2F*)fbTagEffi->Get("eff_udsg_M_ALL");

    //BTag calibration
//     string csv_name_161718[4] = {"DeepCSV_2016LegacySF_V1.csv", "DeepCSV_2016LegacySF_V1.csv", "DeepCSV_106XUL17SF_V2p1.csv", "DeepCSV_106XUL18SF.csv"};
    string csv_name_161718[4] = {"DeepCSV_106XUL16preVFPSF_v1_hzz.csv", "DeepCSV_106XUL16postVFPSF_v2_hzz.csv", "wp_deepCSV_106XUL17_v3_hzz.csv", "wp_deepCSV_106XUL18_v2_hzz.csv"};
    edm::FileInPath btagfileInPath(("Hcc/HccAna/data/"+csv_name_161718[YEAR]).c_str());


bool validate = true; // HARDCODED --> IT COULD BE FALSE!!!
    BTagCalibration calib("DeepCSV", btagfileInPath.fullPath().c_str(), validate);
    reader = new BTagCalibrationReader(BTagEntry::OP_MEDIUM,  // operating point
                                       "central",             // central sys type
                                       {"up", "down"});      // other sys types
   

    reader->load(calib,                // calibration instance
                BTagEntry::FLAV_B,    // btag flavour
                "comb");               // measurement type

    if(year==2018)    {EleBDT_name_161718 = "ElectronMVAEstimatorRun2Fall17IsoV1Values"; BTagCut=0.4184; heepID_name_161718 = "heepElectronID-HEEPV70";}
    if(year==2017)    {EleBDT_name_161718 = "ElectronMVAEstimatorRun2Summer17ULIdIsoValues"; BTagCut=0.4941; heepID_name_161718 = "heepElectronID-HEEPV70";}
    if(year==20165 || year==20160)    {EleBDT_name_161718 = "ElectronMVAEstimatorRun2Summer16ULIdIsoValues"; BTagCut=0.6321; heepID_name_161718 = "heepElectronID-HEEPV70";}*/

}



HccAna::~HccAna()
{
    //destructor --- don't do anything here  
}


// ------------ method called for each event  ------------
void
HccAna::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    using namespace edm;
    using namespace std;
    using namespace pat;
    using namespace trigger;
    using namespace EwkCorrections;

// if(iEvent.id().event() > 709310) 
// 	std::cout<<"PIPPO\t"<<iEvent.id().event()<<"\n";
        
    nEventsTotal += 1.0;

    Run = iEvent.id().run();
    Event = iEvent.id().event();
    LumiSect = iEvent.id().luminosityBlock();

    if (verbose) {
       cout<<"Run: " << Run << ",Event: " << Event << ",LumiSect: "<<LumiSect<<endl;
    }

    // ======= Get Collections ======= //
    if (verbose) {cout<<"getting collections"<<endl;}

    // trigger collection
    edm::Handle<edm::TriggerResults> trigger;
    iEvent.getByToken(triggerSrc_,trigger);
    const edm::TriggerNames trigNames = iEvent.triggerNames(*trigger);

    // trigger Objects
    //edm::Handle<pat::TriggerObjectStandAlone> triggerObjects;
    Handle<std::vector<pat::TriggerObjectStandAlone> > triggerObjects;
    iEvent.getByToken(triggerObjects_, triggerObjects);

    // vertex collection
    edm::Handle<reco::VertexCollection> vertex;
    iEvent.getByToken(vertexSrc_,vertex);

    // electron collection
    edm::Handle<edm::View<pat::Electron> > electrons;
    iEvent.getByToken(elecSrc_,electrons);
    if (verbose) cout<<electrons->size()<<" total electrons in the collection"<<endl;

    // electron before scale/smearing corrections
    edm::Handle<edm::View<pat::Electron> > electronsUnS;
    iEvent.getByToken(elecUnSSrc_,electronsUnS);

    // muon collection
    edm::Handle<edm::View<pat::Muon> > muons;
    iEvent.getByToken(muonSrc_,muons);
    if (verbose) cout<<muons->size()<<" total muons in the collection"<<endl;

    // photon collection 
    /*edm::Handle<edm::View<pat::Photon> > photons;
    iEvent.getByToken(photonSrc_,photons);
    if (verbose) cout<<photons->size()<<" total photons in the collection"<<endl;*/
  
    // met collection 
    edm::Handle<edm::View<pat::MET> > mets;
    iEvent.getByToken(metSrc_,mets);
    
    // Rho Correction
    /*edm::Handle<double> eventRhoMu;
    iEvent.getByToken(muRhoSrc_,eventRhoMu);
    muRho = *eventRhoMu;

    edm::Handle<double> eventRhoE;
    iEvent.getByToken(elRhoSrc_,eventRhoE);
    elRho = *eventRhoE;

    edm::Handle<double> eventRhoSUS;
    iEvent.getByToken(rhoSrcSUS_,eventRhoSUS);
    rhoSUS = *eventRhoSUS;*/

    // Conversions
    edm::Handle< std::vector<reco::Conversion> > theConversions;
    iEvent.getByToken(conversionSrc_, theConversions);
 
    // Beam Spot
    edm::Handle<reco::BeamSpot> beamSpot;
    iEvent.getByToken(beamSpotSrc_,beamSpot);
    const reco::BeamSpot BS = *beamSpot;

    // Particle Flow Cands
    edm::Handle<pat::PackedCandidateCollection> pfCands;
    iEvent.getByToken(pfCandsSrc_,pfCands);

    // FSR Photons
    edm::Handle<edm::View<pat::PFParticle> > photonsForFsr;
    iEvent.getByToken(fsrPhotonsSrc_,photonsForFsr);
  
    // Jets
    //edm::Handle<edm::View<pat::Jet> > jets;
    //iEvent.getByToken(jetSrc_,jets);
		
    // Puppi AK4jets with JESUp
    //edm::Handle<edm::View<pat::Jet> > AK4PuppiJetsJESUp;
    //iEvent.getByToken(AK4PuppiJetJESUpSrc_ ,AK4PuppiJetsJESUp);

    // Puppi AK4jets with JESUp
    //edm::Handle<edm::View<pat::Jet> > AK4PuppiJetsJESDown;
    //iEvent.getByToken(AK4PuppiJetJESDownSrc_ ,AK4PuppiJetsJESDown);

    // Puppi AK4jets with ParticleNet taggers
    edm::Handle<edm::View<pat::Jet> > AK4PuppiJets;
    iEvent.getByToken(AK4PuppiJetSrc_ ,AK4PuppiJets);

    // Puppi AK4jets with ParticleNet taggers
    edm::Handle<edm::View<pat::Jet> > AK4PuppiJetsUncorr;
    iEvent.getByToken(AK4PuppiJetUncorrSrc_ ,AK4PuppiJetsUncorr);
    
    // Puppi AK4jets  Smeared
    edm::Handle<edm::View<pat::Jet> > AK4PuppiJetsSmear;
    iEvent.getByToken(AK4PuppiJetSmearSrc_ ,AK4PuppiJetsSmear);

    // Puppi AK4jets with SmearUp
    edm::Handle<edm::View<pat::Jet> > AK4PuppiJetsSmearUp;
    iEvent.getByToken(AK4PuppiJetSmearUpSrc_ ,AK4PuppiJetsSmearUp);

    // Puppi AK4jets with SmearUp
    edm::Handle<edm::View<pat::Jet> > AK4PuppiJetsSmearDown;
    iEvent.getByToken(AK4PuppiJetSmearDownSrc_ ,AK4PuppiJetsSmearDown);

	// Puppi AK8jets with ParticleNet taggers
	     edm::Handle<edm::View<pat::Jet> > AK8PuppiJets;
	    iEvent.getByToken(AK8PuppiJetSrc_ ,AK8PuppiJets);
	
    //L1 Jets                                       
    //edm::Handle<BXVector<l1t::Jet>> bxvCaloJets;
    //iEvent.getByToken(bxvCaloJetSrc_,bxvCaloJets);
	 
    //L1 Muons
    //edm::Handle<BXVector<l1t::Muon>> bxvCaloMuons;
    //iEvent.getByToken(bxvCaloMuonSrc_,bxvCaloMuons);
		 
    //L1 HT Sum                                       
    //edm::Handle<BXVector<l1t::EtSum>> bxvCaloHT;
    //iEvent.getByToken(bxvCaloHTSrc_,bxvCaloHT);

    //HLT hltAK4PFJetsCorrectedSrc
    /*edm::Handle<edm::View<reco::PFJet>>  hltAK4PFJetsCorrected;
    iEvent.getByToken(hltAK4PFJetsCorrectedSrc_, hltAK4PFJetsCorrected);*/

    /*if (!hltAK4PFJetsCorrected.isValid()) {
        edm::LogWarning("ParticleNetJetTagMonitor") << "Jet collection not valid, will skip the event \n";
        return;
    }*/
		
    //HLT jet for B Tag
    /*edm::Handle<edm::View<reco::PFJet>>  hltjetsForBTag;
    iEvent.getByToken(hltPFJetForBtagSrc_, hltjetsForBTag);

    if (!hltjetsForBTag.isValid()) {
      edm::LogWarning("ParticleNetJetTagMonitor") << "Jet collection not valid, will skip the event \n";
      return;
    }*/
                                
    //edm::Handle<std::vector<reco::PFJet>>  hltjets;
    //iEvent.getByToken(hltPFJetForBtagSrc_, hltjets);

    /*edm::Handle<reco::JetTagCollection> pfJetTagCollectionParticleNetprobc;
    iEvent.getByToken(pfJetTagCollectionParticleNetprobcSrc_, pfJetTagCollectionParticleNetprobc);

    if (!pfJetTagCollectionParticleNetprobc.isValid()) {
      edm::LogWarning("ParticleNetJetTagMonitor") << "HLT Jet tags collection not valid, will skip event \n";
      return;
    }
		
    edm::Handle<reco::JetTagCollection> pfJetTagCollectionParticleNetprobb;
    iEvent.getByToken(pfJetTagCollectionParticleNetprobbSrc_, pfJetTagCollectionParticleNetprobb);

    if (!pfJetTagCollectionParticleNetprobb.isValid()) {
      edm::LogWarning("ParticleNetJetTagMonitor") << "HLT Jet tags collection not valid, will skip event \n";
      return;
    }
		
    edm::Handle<reco::JetTagCollection> pfJetTagCollectionParticleNetprobuds;
    iEvent.getByToken(pfJetTagCollectionParticleNetprobudsSrc_, pfJetTagCollectionParticleNetprobuds);

    if (!pfJetTagCollectionParticleNetprobuds.isValid()) {
      edm::LogWarning("ParticleNetJetTagMonitor") << "HLT Jet tags collection not valid, will skip event \n";
      return;
    }
		
    edm::Handle<reco::JetTagCollection> pfJetTagCollectionParticleNetprobg;
    iEvent.getByToken(pfJetTagCollectionParticleNetprobgSrc_, pfJetTagCollectionParticleNetprobg);

    if (!pfJetTagCollectionParticleNetprobg.isValid()) {
      edm::LogWarning("ParticleNetJetTagMonitor") << "HLT Jet tags collection not valid, will skip event \n";
      return;
    }
		
    edm::Handle<reco::JetTagCollection> pfJetTagCollectionParticleNetprobtauh;
    iEvent.getByToken(pfJetTagCollectionParticleNetprobtauhSrc_, pfJetTagCollectionParticleNetprobtauh);

    if (!pfJetTagCollectionParticleNetprobtauh.isValid()) {
      edm::LogWarning("ParticleNetJetTagMonitor") << "HLT Jet tags collection not valid, will skip event \n";
      return;
    }*/
    /*if (iEvent.getByToken(hltPFJetForBtagSrc_, hltjets)) {
    //get PF jet tags
       edm::Handle<reco::JetTagCollection> pfJetTagCollection;
       bool haveJetTags = false;
       if (iEvent.getByToken(pfJetTagCollectionSrc_, pfJetTagCollection)) {
         haveJetTags = true;
       }
	  }
      cout<<"haveJetTags"<<haveJetTags<<endl;*/
		
    /*if (!jecunc) {



//        edm::ESHandle<JetCorrectorParametersCollection> jetCorrParameterSet;
//        iSetup.get<JetCorrectionsRecord>().get("AK4PFchs", jetCorrParameterSet);


auto const& jetCorrParameterSet = iSetup.getData(mPayloadToken);//"AK4PFchs");
std::vector<JetCorrectorParametersCollection::key_type> keys;
jetCorrParameterSet.validKeys(keys);


//        const JetCorrectorParameters& jetCorrParameters = (*jetCorrParameterSet)["Uncertainty"]; 
        JetCorrectorParameters jetCorrParameters = (jetCorrParameterSet)["Uncertainty"];





        jecunc.reset(new JetCorrectionUncertainty(jetCorrParameters));
    }*/


    // GEN collections
    edm::Handle<reco::GenParticleCollection> prunedgenParticles;
    iEvent.getByToken(prunedgenParticlesSrc_, prunedgenParticles);

    edm::Handle<edm::View<pat::PackedGenParticle> > packedgenParticles;
    iEvent.getByToken(packedgenParticlesSrc_, packedgenParticles);
    
    edm::Handle<edm::View<reco::GenJet> > genJets;
    iEvent.getByToken(genJetsSrc_, genJets);
    
    edm::Handle<GenEventInfoProduct> genEventInfo;
    iEvent.getByToken(generatorSrc_,genEventInfo);
    
    //vector<edm::Handle<LHEEventProduct> > lheInfos;
    //iEvent.getManyByType(lheInfos); // using this method because the label is not always the same (e.g. "source" in the ttH sample)

    edm::Handle<LHEEventProduct> lheInfo;
    iEvent.getByToken(lheInfoSrc_, lheInfo);

    //L1 Trigger
   /* edm::Handle<BXVector<GlobalAlgBlk>> uGtAlgs;
    iEvent.getByToken(algTok_, uGtAlgs);

    if (!uGtAlgs.isValid()) {
        cout << "Cannot find uGT readout record." << endl;
    }*/


//    if (isMC) {    
//        edm::Handle< double > theprefweight;
//            iEvent.getByToken(prefweight_token_, theprefweight ) ;
 //               prefiringWeight =(*theprefweight);
//    }
//    else
        prefiringWeight =1.0;
    
    // ============ Initialize Variables ============= //

    // Event Variables
    if (verbose) {cout<<"clear variables"<<endl;}



    nVtx = -1.0; nInt = -1.0;
    finalState = -1;
    triggersPassed="";
		puN=-1;
    passedTrig=false; passedFullSelection=false; passedZ4lSelection=false; passedQCDcut=false; passedZqqSelection=false;
    Trigger_l1name.clear();
    Trigger_l1decision.clear();
    Trigger_hltname.clear();
    Trigger_hltdecision.clear();

    // Event Weights
    genWeight=1.0; pileupWeight=1.0; pileupWeightUp=1.0; pileupWeightDn=1.0; dataMCWeight=1.0; eventWeight=1.0;
    k_qqZZ_qcd_dPhi = 1.0; k_qqZZ_qcd_M = 1.0; k_qqZZ_qcd_Pt = 1.0; k_qqZZ_ewk = 1.0;

    qcdWeights.clear(); nnloWeights.clear(); pdfWeights.clear();
    pdfRMSup=1.0; pdfRMSdown=1.0; pdfENVup=1.0; pdfENVdown=1.0;

    //lepton variables
	
	//ALLlep_pt.clear(); ALLlep_eta.clear(); ALLlep_phi.clear(); ALLlep_mass.clear(); ALLlep_id.clear();
	Ele_pt.clear(); Ele_eta.clear(); Ele_phi.clear(); Ele_mass.clear(); Ele_dxy.clear(); Ele_dz.clear(); Ele_id.clear(); Ele_hcalIso.clear(); Ele_ecalIso.clear(); Ele_trackIso.clear(); Ele_isEB.clear(); Ele_IsoCal.clear(); /*Ele_PF_Iso_R04.clear();*/ Ele_isPassID.clear(); Ele_03_Neutral.clear();
    Muon_pt.clear(); Muon_eta.clear(); Muon_phi.clear(); Muon_mass.clear(); Muon_dxy.clear(); Muon_dz.clear(); Muon_id.clear(); Muon_isPF.clear(); Muon_PF_Iso_R04.clear(); Muon_PassLooseID.clear();
    AK4lep_pt.clear(); AK4lep_eta.clear(); AK4lep_phi.clear(); AK4lep_mass.clear(); AK4lep_id.clear();
		
    //hlt Jets for b tag
    hltjetForBTag_pt.clear();
    hltjetForBTag_eta.clear();
    hltjetForBTag_phi.clear();
    hltjetForBTag_mass.clear();
    hltParticleNetONNXJetTags_probb.clear(); hltParticleNetONNXJetTags_probc.clear(); hltParticleNetONNXJetTags_probuds.clear(); hltParticleNetONNXJetTags_probg.clear(); hltParticleNetONNXJetTags_probtauh.clear();
	
    //hltAK4PFJetsCorrected    
		
    hltAK4PFJetsCorrected_pt.clear();
    hltAK4PFJetsCorrected_eta.clear();
    hltAK4PFJetsCorrected_phi.clear();
    hltAK4PFJetsCorrected_mass.clear();

    //HLT jets for turn on curves and scale factors
    HLTJet80_pt.clear();
    HLTJet80_eta.clear();
    HLTJet80_phi.clear();

    HLTJet60_pt.clear();
    HLTJet60_eta.clear();
    HLTJet60_phi.clear();

    HLTJet80_MatchedCalo_pt.clear();
    HLTJet80_MatchedCalo_eta.clear();
    HLTJet80_MatchedCalo_phi.clear();

    HLTJet60_MatchedCalo_pt.clear();
    HLTJet60_MatchedCalo_eta.clear();
    HLTJet60_MatchedCalo_phi.clear();

    HLTAK4PFJetLoose_pt.clear();
    HLTAK4PFJetLoose_eta.clear();
    HLTAK4PFJetLoose_phi.clear();

    HLTAK4PFJetTight_pt.clear();
    HLTAK4PFJetTight_eta.clear();
    HLTAK4PFJetTight_phi.clear();

    HLTAK4PFJet_pt.clear();
    HLTAK4PFJet_eta.clear();
    HLTAK4PFJet_phi.clear();
    // Puppi AK4jets with ParticleNet taggers
    AK4PuppiJets_pt.clear();
    AK4PuppiJetsJESUp_pt.clear();
    AK4PuppiJetsJESDown_pt.clear();
    AK4PuppiJetsSmear_pt.clear();
    AK4PuppiJetsSmearUp_pt.clear();
    AK4PuppiJetsSmearDown_pt.clear();
    AK4PuppiJets_uncorrpt.clear();
    AK4PuppiJets_eta.clear();
    AK4PuppiJets_phi.clear();
    AK4PuppiJets_mass.clear();
    AK4PuppiJets_hadrFlav.clear();
    AK4PuppiJets_partFlav.clear();
    AK4PuppiJets_partFlav.clear();

    jet_pfParticleNetAK4JetTags_probb.clear(); jet_pfParticleNetAK4JetTags_probc.clear(); jet_pfParticleNetAK4JetTags_probuds.clear(); jet_pfParticleNetAK4JetTags_probg.clear(); jet_pfParticleNetAK4JetTags_probtauh.clear();
    jet_pfParticleNetAK4JetTags_CvsB.clear(); jet_pfParticleNetAK4JetTags_CvsL.clear(); jet_pfParticleNetAK4JetTags_CvsAll.clear(); jet_pfParticleNetAK4JetTags_BvsC.clear(); jet_pfParticleNetAK4JetTags_BvsL.clear(); jet_pfParticleNetAK4JetTags_BvsAll.clear();  jet_pfParticleNetAK4JetTags_QvsG.clear();
    jet_PNetRegPtRawCorr.clear(); jet_PNetRegPtRawCorrNeutrino.clear(); jet_PNetRegPtRawRes.clear();





    jet_pfDeepJetAK4JetTags_probb.clear(); jet_pfDeepJetAK4JetTags_probbb.clear(); jet_pfDeepJetAK4JetTags_problepb.clear(); jet_pfDeepJetAK4JetTags_probc.clear(); jet_pfDeepJetAK4JetTags_probuds.clear(); jet_pfDeepJetAK4JetTags_probg.clear();

    jet_pfDeepCSVAK4JetTags_probb.clear(); jet_pfDeepCSVAK4JetTags_probbb.clear(); jet_pfDeepCSVAK4JetTags_probc.clear(); jet_pfDeepCSVAK4JetTags_probudsg.clear();
    
    // Puppi AK8jets with ParticleNet and DeepDoubleX taggers
    leadingAK8_pt_idx = -1;
    subleadingAK8_pt_idx = -1;

    AK8PuppiJets_pt.clear();
    AK8PuppiJets_eta.clear();
    AK8PuppiJets_phi.clear();
    AK8PuppiJets_mass.clear();
    AK8PuppiJets_softdropmass.clear();

    jet_pfParticleNetJetTags_probZbb.clear(); jet_pfParticleNetJetTags_probZcc.clear(); jet_pfParticleNetJetTags_probZqq.clear(); jet_pfParticleNetJetTags_probQCDbb .clear();  jet_pfParticleNetJetTags_probQCDcc.clear(); jet_pfParticleNetJetTags_probQCDb.clear(); jet_pfParticleNetJetTags_probQCDc.clear(); jet_pfParticleNetJetTags_probQCDothers.clear(); jet_pfParticleNetJetTags_probHbb.clear(); jet_pfParticleNetJetTags_probHcc.clear(); jet_pfParticleNetJetTags_probHqqqq.clear(); 
    jet_pfMassDecorrelatedParticleNetJetTags_probXbb.clear(); jet_pfMassDecorrelatedParticleNetJetTags_probXcc.clear(); jet_pfMassDecorrelatedParticleNetJetTags_probXqq.clear(); jet_pfMassDecorrelatedParticleNetJetTags_probQCDbb.clear(); jet_pfMassDecorrelatedParticleNetJetTags_probQCDcc.clear(); jet_pfMassDecorrelatedParticleNetJetTags_probQCDb.clear(); jet_pfMassDecorrelatedParticleNetJetTags_probQCDc.clear(); jet_pfMassDecorrelatedParticleNetJetTags_probQCDothers.clear();
    jet_pfMassIndependentDeepDoubleBvLV2JetTags_probHbb.clear(); jet_pfMassIndependentDeepDoubleCvLV2JetTags_probHcc.clear(); jet_pfMassIndependentDeepDoubleCvBV2JetTags_probHcc.clear();
    
    // MET
    met=-1.0;  met_pt=-1.0; met_phi=9999.0;
    //met_jesup=-1.0; met_phi_jesup=9999.0; met_jesdn=-1.0; met_phi_jesdn=9999.0; 
    //met_uncenup=-1.0; met_phi_uncenup=9999.0; met_uncendn=-1.0; met_phi_uncendn=9999.0; 

    L1jet_pt.clear(); L1jet_eta.clear(); L1jet_phi.clear(); L1jet_mass.clear();
    L1muon_pt.clear(); L1muon_eta.clear(); L1muon_phi.clear(); L1muon_mass.clear(); L1muon_qual.clear();

    // -------------------------
    // GEN level information
    // ------------------------- 

    //Event variables
    GENfinalState=-1;

    // Jets
    GENjet_pt.clear(); GENjet_eta.clear(); GENjet_phi.clear(); GENjet_mass.clear(); 
    GENnjets_pt30_eta4p7=0;
    GENnjets_pt30_eta2p5=0;
    GENpt_leadingjet_pt30_eta4p7=-1.0; GENabsrapidity_leadingjet_pt30_eta4p7=-1.0; GENabsdeltarapidity_hleadingjet_pt30_eta4p7=-1.0;
    GENpt_leadingjet_pt30_eta2p5=-1.0; 
    lheNb=0; lheNj=0; nGenStatus2bHad=0;

    //quarks
    quark_pt.clear(); quark_eta.clear(); quark_phi.clear(); quark_flavour.clear(); quark_VBF.clear();
    quark_pt_float.clear(); quark_eta_float.clear(); quark_phi_float.clear();
    Z_pt = -1000; Z_eta = -1000; Z_phi = -1000; Z_mass = -1000;
    //

    if (verbose) {cout<<"clear other variables"<<endl; }
    // Resolution
    //massErrorUCSD=-1.0; massErrorUCSDCorr=-1.0; massErrorUF=-1.0; massErrorUFCorr=-1.0; massErrorUFADCorr=-1.0;

    // Event Category
    EventCat=-1;

    hltjetForBTag_pt_float.clear(); hltjetForBTag_eta_float.clear(); hltjetForBTag_phi_float.clear(); hltjetForBTag_mass_float.clear();
    hltAK4PFJetsCorrected_pt_float.clear(); hltAK4PFJetsCorrected_eta_float.clear(); hltAK4PFJetsCorrected_phi_float.clear(); hltAK4PFJetsCorrected_mass_float.clear();

    L1jet_pt_float.clear(); L1jet_eta_float.clear(); L1jet_phi_float.clear(); L1jet_mass_float.clear();
    L1muon_pt_float.clear(); L1muon_eta_float.clear(); L1muon_phi_float.clear(); L1muon_mass_float.clear();

    AK4PuppiJets_pt_float.clear(); 
    AK4PuppiJetsJESUp_pt_float.clear();AK4PuppiJetsJESDown_pt_float.clear(); 
    AK4PuppiJetsSmear_pt_float.clear(); AK4PuppiJetsSmearUp_pt_float.clear();AK4PuppiJetsSmearDown_pt_float.clear();  AK4PuppiJets_uncorrpt_float.clear(); AK4PuppiJets_eta_float.clear(); AK4PuppiJets_phi_float.clear(); AK4PuppiJets_mass_float.clear();
	AK8PuppiJets_pt_float.clear(); AK8PuppiJets_eta_float.clear(); AK8PuppiJets_phi_float.clear(); AK8PuppiJets_mass_float.clear();
    AK4PuppiJets_qgl.clear();

    HLTJet80_pt_float.clear();  HLTJet80_eta_float.clear(); HLTJet80_phi_float.clear();
    HLTJet60_pt_float.clear();  HLTJet60_eta_float.clear(); HLTJet60_phi_float.clear();
    HLTJet80_MatchedCalo_pt_float.clear(); HLTJet80_MatchedCalo_eta_float.clear(); HLTJet80_MatchedCalo_phi_float.clear(); 
    HLTJet60_MatchedCalo_pt_float.clear(); HLTJet60_MatchedCalo_eta_float.clear(); HLTJet60_MatchedCalo_phi_float.clear(); 
    HLTAK4PFJetLoose_pt_float.clear(); HLTAK4PFJetLoose_eta_float.clear(); HLTAK4PFJetLoose_phi_float.clear();
    HLTAK4PFJetTight_pt_float.clear(); HLTAK4PFJetTight_eta_float.clear(); HLTAK4PFJetTight_phi_float.clear();
    HLTAK4PFJet_pt_float.clear(); HLTAK4PFJet_eta_float.clear(); HLTAK4PFJet_phi_float.clear();

    // ====================== Do Analysis ======================== //
// if(iEvent.id().event() > 709310) 
// 	std::cout<<"PIPPO\tdopo inizializzazione\n";
		//cout<<"aaa"<<endl;
    std::map<int, TLorentzVector> fsrmap;
    vector<reco::Candidate*> selectedLeptons;
    std::map<unsigned int, TLorentzVector> selectedFsrMap;

    fsrmap.clear(); selectedFsrMap.clear(); selectedLeptons.clear();

    if (verbose) cout<<"start pileup reweighting"<<endl;

    // PU information
    if(isMC && reweightForPU) {        
       edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
        iEvent.getByToken(pileupSrc_, PupInfo);
				puN = PupInfo->begin()->getTrueNumInteractions();      

        if (verbose) cout<<"got pileup info"<<endl;

        std::vector<PileupSummaryInfo>::const_iterator PVI;      
        int npv = -1;
        for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
            int BX = PVI->getBunchCrossing();
            if(BX == 0) { npv = PVI->getTrueNumInteractions(); continue;}
        }        
        if (verbose) cout<<"N true interations = "<<npv<<endl;
        nInt = npv;
        //pileupWeight = pileUp.getPUWeight(npv,PUVersion);
        pileupWeight = pileUp.getPUWeight(h_pileup,npv);
//std::cout<<pileupWeight<<"\t"<<npv<<std::endl;
        pileupWeightUp = pileUp.getPUWeight(h_pileupUp,npv);
        pileupWeightDn = pileUp.getPUWeight(h_pileupDn,npv);
        if (verbose) cout<<"pileup weight = "<<pileupWeight<<", filling histograms"<<endl;
        histContainer_["NINTERACT"]->Fill(npv);
        histContainer_["NINTERACT_RW"]->Fill(npv,pileupWeight);
    } else { pileupWeight = 1.0;}   

    if (verbose) {cout<<"finished pileup reweighting"<<endl; }
    
    if(isMC) {
        float tmpWeight = genEventInfo->weight();
        genWeight = (tmpWeight > 0 ? 1.0 : -1.0);
        if (verbose) {cout<<"tmpWeight: "<<tmpWeight<<"; genWeight: "<<genWeight<<endl;}        
        double rms = 0.0;

        //std::cout<<"tmpWeight: "<<tmpWeight<<std::endl;

        if(lheInfo.isValid()){
            
            for(unsigned int i = 0; i < lheInfo->weights().size(); i++) {

                tmpWeight = genEventInfo->weight();
                tmpWeight *= lheInfo->weights()[i].wgt/lheInfo->originalXWGTUP();
                pdfWeights.push_back(tmpWeight);

                if (i<=8 or int(i)>=posNNPDF) {
                    tmpWeight = genEventInfo->weight();
                    tmpWeight *= lheInfo->weights()[i].wgt/lheInfo->originalXWGTUP();
                    if (int(i)<posNNPDF) {qcdWeights.push_back(tmpWeight);}
                }
                else {
                    tmpWeight = lheInfo->weights()[i].wgt;
                    tmpWeight /= lheInfo->originalXWGTUP();
                    //if (i==9) genWeight = tmpWeight;
                    if (int(i)<posNNPDF) {nnloWeights.push_back(tmpWeight);}
                }
                // NNPDF30 variations
                if (int(i)>=posNNPDF && int(i)<=(posNNPDF+100)) {
                    rms += tmpWeight*tmpWeight;
                    if (tmpWeight>pdfENVup) pdfENVup=tmpWeight;
                    if (tmpWeight<pdfENVdown) pdfENVdown=tmpWeight;
                }
            }
            pdfRMSup=sqrt(rms/100.0); pdfRMSdown=1.0/pdfRMSup;
            if (verbose) cout<<"pdfRMSup "<<pdfRMSup<<" pdfRMSdown "<<pdfRMSdown<<endl;
        
            const lhef::HEPEUP& lheEvent = lheInfo->hepeup();
            std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;
            for ( size_t idxParticle = 0; idxParticle < lheParticles.size(); ++idxParticle ) {
                int id = std::abs(lheEvent.IDUP[idxParticle]);
                int status = lheEvent.ISTUP[idxParticle];
                if ( status == 1 && id==5 ) { 
                    lheNb += 1;
                }
                if ( status == 1 && ((id >= 1 && id <= 6) || id == 21) ) { 
                    lheNj += 1;
                }
            }
        
        }
        
        if (verbose) cout<<"setting gen variables"<<endl;       
        setGENVariables(prunedgenParticles,packedgenParticles,genJets); 
        if (verbose) { cout<<"finshed setting gen variables"<<endl;  }


    } //end if isMC
    sumWeightsTotal += genWeight;
    sumWeightsTotalPU += pileupWeight*genWeight;

    eventWeight = pileupWeight*genWeight;
    
    // Fill L1 seeds and decisions
    /*gtUtil_->retrieveL1(iEvent, iSetup, algInputTag_);
    const vector<pair<string, bool> > finalDecisions = gtUtil_->decisionsFinal();
    for (size_t i_l1t = 0; i_l1t < finalDecisions.size(); i_l1t++){
        string l1tName = (finalDecisions.at(i_l1t)).first;
        if( l1tName.find("SingleJet") != string::npos || l1tName.find("DoubleJet") != string::npos || l1tName.find("TripleJet") != string::npos || l1tName.find("QuadJet") != string::npos ||  l1tName.find("HTT")!= string::npos ){
            //cout << "L1: " << l1tName << " | decision: " << finalDecisions.at(i_l1t).second << endl;
            Trigger_l1name.push_back( l1tName );
            Trigger_l1decision.push_back( finalDecisions.at(i_l1t).second );
        }
      }*/
    //bool trigger_PFJet80=false;
    unsigned int _tSize = trigger->size();
    // create a string with all passing trigger names
    for (unsigned int i=0; i<_tSize; ++i) {
        std::string triggerName = trigNames.triggerName(i);
        if (strstr(triggerName.c_str(),"_step")) continue;
        if (strstr(triggerName.c_str(),"MC_")) continue;
        if (strstr(triggerName.c_str(),"AlCa_")) continue;
        if (strstr(triggerName.c_str(),"DST_")) continue;
        if (strstr(triggerName.c_str(),"HLT_HI")) continue;
        if (strstr(triggerName.c_str(),"HLT_Physics")) continue;
        if (strstr(triggerName.c_str(),"HLT_Random")) continue;
        if (strstr(triggerName.c_str(),"HLT_ZeroBias")) continue;
        if (strstr(triggerName.c_str(),"HLT_IsoTrack")) continue;
        if (strstr(triggerName.c_str(),"Hcal")) continue;
        if (strstr(triggerName.c_str(),"Ecal")) continue;
        if (trigger->accept(i)) triggersPassed += triggerName+" ";
			
				
        //if(triggerName.find("HLT_QuadPFJet70_50_45_35_PFBTagParticleNet_2BTagSum0p65") != string::npos || triggerName.find("HLT_PFJet500") != string::npos ){
        //if(triggerName.find("HLT_QuadPFJet") != string::npos || triggerName.find("HLT_PFJet") != string::npos || triggerName.find("HLT_DiPFJetAve") != string::npos || triggerName.find("HLT_AK8PFJet") != string::npos ) {
        Trigger_hltname.push_back(triggerName);
        Trigger_hltdecision.push_back(trigger->accept(i));
        //if(triggerName.find("HLT_PFJet80_v") != string::npos && trigger->accept(i)==1){
        //  trigger_PFJet80=true;
        //}
    //}
    }
    if (firstEntry) cout<<"triggersPassed: "<<triggersPassed<<endl;
    //cout<<"triggersPassed: "<<triggersPassed<<endl;
    firstEntry = false;
    // check if any of the triggers in the user list have passed
    bool passedSingleEl=false;
    bool passedSingleMu=false;
    bool passedAnyOther=false;
    for (unsigned int i=0; i<triggerList.size(); ++i) {
        if (strstr(triggersPassed.c_str(),triggerList.at(i).c_str())) {
            passedTrig=true;
            if (!isMC) {
								//cout<<"bbb"<<endl;
                if (strstr(triggerList.at(i).c_str(),"_WP")) passedSingleEl=true;
                if (strstr(triggerList.at(i).c_str(),"HLT_Iso")) passedSingleMu=true;
                if (strstr(triggerList.at(i).c_str(),"CaloIdL")) passedAnyOther=true;
                if (strstr(triggerList.at(i).c_str(),"TrkIsoVVL")) passedAnyOther=true;
                if (strstr(triggerList.at(i).c_str(),"Triple")) passedAnyOther=true;
            }
        }
    }
    
    bool passedOnlySingle=((passedSingleEl && !passedAnyOther) || (passedSingleMu && !passedSingleEl && !passedAnyOther));
    bool trigConditionData = ( passedTrig && (!checkOnlySingle || (checkOnlySingle && passedOnlySingle)) );
if(trigConditionData && verbose)
	std::cout<<""<<std::endl;

//    bool trigConditionData = true;
        
    if (verbose) cout<<"checking PV"<<endl;       
    const reco::Vertex *PV = 0;
    int theVertex = -1;
    for (unsigned int i=0; i<vertex->size(); i++) {
        PV = &(vertex->at(i));        
        if (verbose) std::cout<<"isFake: "<<PV->isFake()<<" chi2 "<<PV->chi2()<<" ndof "<<PV->ndof()<<" rho "<<PV->position().Rho()<<" Z "<<PV->position().Z()<<endl; 
        //if (PV->chi2()==0 && PV->ndof()==0) continue;
        if (PV->isFake()) continue;
        if (PV->ndof()<=4 || PV->position().Rho()>2.0 || fabs(PV->position().Z())>24.0) continue;
        theVertex=(int)i; break;
    }        

    if (verbose) std::cout<<"vtx: "<<theVertex<<" trigConditionData "<<trigConditionData<<" passedTrig "<<passedTrig<<std::endl;
 
    //if(theVertex >= 0 && (isMC || (!isMC && trigConditionData)) )  {
    if(theVertex >= 0 && (isMC || (!isMC )) )  {

        if (verbose) cout<<"good PV "<<theVertex<<endl; 
        
        PV_x =  PV->position().X();
        PV_y =  PV->position().Y();
        PV_z =  PV->position().Z();

        BS_x =  BS.position().X();
        BS_y =  BS.position().Y();
        BS_z =  BS.position().Z();
        BS_xErr =  BS.x0Error();
        BS_yErr =  BS.y0Error();
        BS_zErr =  BS.z0Error();

        BeamWidth_x = BS.BeamWidthX();
        BeamWidth_y = BS.BeamWidthY();
        BeamWidth_xErr = BS.BeamWidthXError();
        BeamWidth_yErr = BS.BeamWidthYError();
    
        //N Vertex 
        if (verbose) {cout<<"fill nvtx histogram"<<endl;}
        nVtx = vertex->size();
        histContainer_["NVTX"]->Fill(nVtx);
        histContainer_["NVTX_RW"]->Fill(nVtx,pileupWeight);
// if(iEvent.id().event() > 709310){
// 	std::cout<<"PIPPO\tdopo vertex info\n";
// }
        //MET
        if (verbose) {cout<<"get met value"<<endl;}
        if (!mets->empty()) {
            met = (*mets)[0].et();
            met_pt = (*mets)[0].pt();
            met_phi = (*mets)[0].phi();
            /*met_jesup = (*mets)[0].shiftedPt(pat::MET::JetEnUp);
            met_phi_jesup = (*mets)[0].shiftedPhi(pat::MET::JetEnUp);
            met_jesdn = (*mets)[0].shiftedPt(pat::MET::JetEnDown);
            met_phi_jesdn = (*mets)[0].shiftedPhi(pat::MET::JetEnDown);
            met_uncenup = (*mets)[0].shiftedPt(pat::MET::UnclusteredEnUp);
            met_phi_uncenup = (*mets)[0].shiftedPhi(pat::MET::UnclusteredEnUp);
            met_uncendn = (*mets)[0].shiftedPt(pat::MET::UnclusteredEnDown);
            met_phi_uncendn = (*mets)[0].shiftedPhi(pat::MET::UnclusteredEnDown);*/        
        }

        if (verbose) cout<<"start lepton analysis"<<endl;           
        vector<pat::Electron> AllElectrons; 
        vector<pat::Electron> AllElectronsUnS;////uncorrected electron 
        vector<pat::Muon> AllMuons; 
       // vector<pat::Tau> AllTaus; 
       // vector<pat::Photon> AllPhotons;
        AllElectrons = helper.goodLooseElectrons2012(electrons,_elecPtCut);
       // AllElectronsUnS = helper.goodLooseElectrons2012(electrons,electronsUnS,_elecPtCut);
        AllMuons = helper.goodLooseMuons2012(muons,_muPtCut);
       // AllTaus = helper.goodLooseTaus2015(taus,_tauPtCut);
       // AllPhotons = helper.goodLoosePhotons2015(photons,_phoPtCut);

        /*helper.cleanOverlappingLeptons(AllMuons,AllElectrons,PV);
        helper.cleanOverlappingLeptons(AllMuons,AllElectronsUnS,PV);
        recoMuons = helper.goodMuons2015_noIso_noPf(AllMuons,_muPtCut,PV);
        recoElectrons = helper.goodElectrons2015_noIso_noBdt(AllElectrons,_elecPtCut,elecID,PV,iEvent,sip3dCut, true);
        recoElectronsUnS = helper.goodElectrons2015_noIso_noBdt(AllElectronsUnS,_elecPtCut,elecID,PV,iEvent,sip3dCut, false);
        helper.cleanOverlappingTaus(recoMuons,recoElectrons,AllTaus,isoCutMu,isoCutEl,muRho,elRho);
        recoTaus = helper.goodTaus2015(AllTaus,_tauPtCut);
        recoPhotons = helper.goodPhotons2015(AllPhotons,_phoPtCut,year);*/
       
        if (verbose) cout<<"before vector assign"<<std::endl;
        std::cout<<"event: "<<iEvent.id().event()<<std::endl;
	    setTreeVariables(iEvent, iSetup,AK4PuppiJets, AK4PuppiJetsUncorr, uncAK4,AK4PuppiJetsSmear, AK4PuppiJetsSmearUp,AK4PuppiJetsSmearDown, AK8PuppiJets,  AllMuons, AllElectrons, PV);
		//setTreeVariables(iEvent, iSetup, goodJets, selectedMergedJets, AK4PuppiJets, AK8PuppiJets, bxvCaloJets, bxvCaloMuons, bxvCaloHT, AllMuons, AllElectrons, PV);		
        //setTreeVariables(iEvent, iSetup, goodJets, goodJetQGTagger,goodJetaxis2, goodJetptD, goodJetmult, selectedMergedJets, hltjetsForBTag,  hltAK4PFJetsCorrected, pfJetTagCollectionParticleNetprobc , pfJetTagCollectionParticleNetprobb , pfJetTagCollectionParticleNetprobuds , pfJetTagCollectionParticleNetprobg ,pfJetTagCollectionParticleNetprobtauh ,  bxvCaloJets, bxvCaloMuons, bxvCaloHT, AllMuons, AllElectrons);

      	if (verbose) cout<<"finshed setting tree variables"<<endl;


        //HLT jets for turn on curves and scale factors
        vector<pat::TriggerObjectStandAlone> TriggerObj_PFJet60, TriggerObj_PFJet80, TriggerObj_AK4PFJetsLoose, TriggerObj_AK4PFJetsTight, TriggerObj_AK4PFJets, TriggerObj_PFJet80MatchedCalo, TriggerObj_PFJet60MatchedCalo;
        const TriggerNames &triggerNames = iEvent.triggerNames( *trigger );

        for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
          //std::cout<<"cycle on trigger objects"<<std::endl;
          obj.unpackPathNames(triggerNames);

          if(obj.collection()=="hltAK4PFJetsTightIDCorrected::HLT"){
            TriggerObj_AK4PFJetsTight.push_back(obj);
          }

          if(obj.collection()=="hltAK4PFJetsLooseIDCorrected::HLT"){
            TriggerObj_AK4PFJetsLoose.push_back(obj);
          }

          if(obj.collection()=="hltAK4PFJetsCorrected::HLT"){
            TriggerObj_AK4PFJets.push_back(obj);
          }

          if(obj.collection()=="hltPFJetsCorrectedMatchedToCaloJets50::HLT"){
            TriggerObj_PFJet80MatchedCalo.push_back(obj);
          }

          if(obj.collection()=="hltPFJetsCorrectedMatchedToCaloJets40::HLT"){
            TriggerObj_PFJet60MatchedCalo.push_back(obj);
          }
          for (unsigned h = 0; h < obj.filterLabels().size(); ++h) {
            //if(trigger_PFJet80==true){
            //  cout<<"Trig objs: "<< obj.filterLabels()[h]<<endl;
            //}
            //std::cout<<"cycle on filterLabels"<<std::endl;
            if(obj.filterLabels()[h]=="hltSinglePFJet80"){
               //cout<<"BBBB"<<endl; 
               TriggerObj_PFJet80.push_back(obj);
            }

            if(obj.filterLabels()[h]=="hltSinglePFJet60"){
               //cout<<"BBBB"<<endl; 
               TriggerObj_PFJet60.push_back(obj);
            }
            /*if(obj.filterLabels()[h]=="hltPFJetsCorrectedMatchedToCaloJets50"){
                std::cout<<"AAAAA"<<std::endl;
                TriggerObj_AK4PFJetsLoose.push_back(obj);
            }

            if(obj.filterLabels()[h]=="hltAK4PFJetsTightIDCorrected"){
                TriggerObj_AK4PFJetsTight.push_back(obj);
            }*/
          }

        }  
  
        for(uint t=0; t<TriggerObj_PFJet80.size(); t++){
          HLTJet80_pt.push_back(TriggerObj_PFJet80.at(t).pt());
          HLTJet80_eta.push_back(TriggerObj_PFJet80.at(t).eta());
          HLTJet80_phi.push_back(TriggerObj_PFJet80.at(t).phi());
          //cout<<" hlt jet pt (80): "<<TriggerObj_PFJet80.at(t).pt()<<endl;
        }

        for(uint t=0; t<TriggerObj_PFJet60.size(); t++){
          HLTJet60_pt.push_back(TriggerObj_PFJet60.at(t).pt());
          HLTJet60_eta.push_back(TriggerObj_PFJet60.at(t).eta());
          HLTJet60_phi.push_back(TriggerObj_PFJet60.at(t).phi());
          //cout<<" hlt jet pt (80): "<<TriggerObj_PFJet80.at(t).pt()<<endl;
        }

        for(uint t=0; t<TriggerObj_PFJet80MatchedCalo.size(); t++){
          HLTJet80_MatchedCalo_pt.push_back(TriggerObj_PFJet80MatchedCalo.at(t).pt());
          HLTJet80_MatchedCalo_eta.push_back(TriggerObj_PFJet80MatchedCalo.at(t).eta());
          HLTJet80_MatchedCalo_phi.push_back(TriggerObj_PFJet80MatchedCalo.at(t).phi());
          //cout<<" hlt jet pt (80): "<<TriggerObj_PFJet80.at(t).pt()<<endl;
        }

        for(uint t=0; t<TriggerObj_PFJet60MatchedCalo.size(); t++){
          HLTJet60_MatchedCalo_pt.push_back(TriggerObj_PFJet60MatchedCalo.at(t).pt());
          HLTJet60_MatchedCalo_eta.push_back(TriggerObj_PFJet60MatchedCalo.at(t).eta());
          HLTJet60_MatchedCalo_phi.push_back(TriggerObj_PFJet60MatchedCalo.at(t).phi());
          //cout<<" hlt jet pt (80): "<<TriggerObj_PFJet80.at(t).pt()<<endl;
        }

        for(uint t=0; t<TriggerObj_AK4PFJetsLoose.size(); t++){
          HLTAK4PFJetLoose_pt.push_back(TriggerObj_AK4PFJetsLoose.at(t).pt());
          HLTAK4PFJetLoose_eta.push_back(TriggerObj_AK4PFJetsLoose.at(t).eta());
          HLTAK4PFJetLoose_phi.push_back(TriggerObj_AK4PFJetsLoose.at(t).phi());
          //cout<<" hlt jet pt: "<<TriggerObj_AK4PFJetsLoose.at(t).pt()<<endl;
        }
        
        for(uint t=0; t<TriggerObj_AK4PFJetsTight.size(); t++){
          HLTAK4PFJetTight_pt.push_back(TriggerObj_AK4PFJetsTight.at(t).pt());
          HLTAK4PFJetTight_eta.push_back(TriggerObj_AK4PFJetsTight.at(t).eta());
          HLTAK4PFJetTight_phi.push_back(TriggerObj_AK4PFJetsTight.at(t).phi());
          //cout<<" hlt jet pt: "<<TriggerObj_AK4PFJetsTight.at(t).pt()<<endl;
        }
        for(uint t=0; t<TriggerObj_AK4PFJets.size(); t++){
          HLTAK4PFJet_pt.push_back(TriggerObj_AK4PFJets.at(t).pt());
          HLTAK4PFJet_eta.push_back(TriggerObj_AK4PFJets.at(t).eta());
          HLTAK4PFJet_phi.push_back(TriggerObj_AK4PFJets.at(t).phi());
          //cout<<" hlt jet pt: "<<TriggerObj_AK4PFJetsTight.at(t).pt()<<endl;
        }

        quark_pt_float.assign(quark_pt.begin(),quark_pt.end());
        quark_eta_float.assign(quark_eta.begin(),quark_eta.end());
        quark_phi_float.assign(quark_phi.begin(),quark_phi.end());
	
        L1jet_pt_float.assign(L1jet_pt.begin(),L1jet_pt.end());
        L1jet_eta_float.assign(L1jet_eta.begin(),L1jet_eta.end());
        L1jet_phi_float.assign(L1jet_phi.begin(),L1jet_phi.end());
        L1jet_mass_float.assign(L1jet_mass.begin(),L1jet_mass.end());

        L1muon_pt_float.assign(L1muon_pt.begin(),L1muon_pt.end());
        L1muon_eta_float.assign(L1muon_eta.begin(),L1muon_eta.end());
        L1muon_phi_float.assign(L1muon_phi.begin(),L1muon_phi.end());
      	L1muon_mass_float.assign(L1muon_mass.begin(),L1muon_mass.end());

        hltjetForBTag_pt_float.assign(hltjetForBTag_pt.begin(), hltjetForBTag_pt.end());
        hltjetForBTag_eta_float.assign(hltjetForBTag_eta.begin(), hltjetForBTag_eta.end());
        hltjetForBTag_phi_float.assign(hltjetForBTag_phi.begin(), hltjetForBTag_phi.end());
        hltjetForBTag_mass_float.assign(hltjetForBTag_mass.begin(), hltjetForBTag_mass.end());

				
        hltAK4PFJetsCorrected_pt_float.assign(hltAK4PFJetsCorrected_pt.begin(), hltAK4PFJetsCorrected_pt.end());
        hltAK4PFJetsCorrected_eta_float.assign(hltAK4PFJetsCorrected_eta.begin(), hltAK4PFJetsCorrected_eta.end());
        hltAK4PFJetsCorrected_phi_float.assign(hltAK4PFJetsCorrected_phi.begin(), hltAK4PFJetsCorrected_phi.end());
        hltAK4PFJetsCorrected_mass_float.assign(hltAK4PFJetsCorrected_mass.begin(), hltAK4PFJetsCorrected_mass.end());

        AK4PuppiJets_pt_float.assign(AK4PuppiJets_pt.begin(), AK4PuppiJets_pt.end()); 
        AK4PuppiJetsJESUp_pt_float.assign(AK4PuppiJetsJESUp_pt.begin(), AK4PuppiJetsJESUp_pt.end()); 
        AK4PuppiJetsJESDown_pt_float.assign(AK4PuppiJetsJESDown_pt.begin(), AK4PuppiJetsJESDown_pt.end()); 
        AK4PuppiJetsSmear_pt_float.assign(AK4PuppiJetsSmear_pt.begin(), AK4PuppiJetsSmear_pt.end()); 
        AK4PuppiJetsSmearUp_pt_float.assign(AK4PuppiJetsSmearUp_pt.begin(), AK4PuppiJetsSmearUp_pt.end()); 
        AK4PuppiJetsSmearDown_pt_float.assign(AK4PuppiJetsSmearDown_pt.begin(), AK4PuppiJetsSmearDown_pt.end()); 
        AK4PuppiJets_uncorrpt_float.assign(AK4PuppiJets_uncorrpt.begin(), AK4PuppiJets_uncorrpt.end()); 
        AK4PuppiJets_eta_float.assign(AK4PuppiJets_eta.begin(), AK4PuppiJets_eta.end()); 
        AK4PuppiJets_phi_float.assign(AK4PuppiJets_phi.begin(), AK4PuppiJets_phi.end()); 
        AK4PuppiJets_mass_float.assign(AK4PuppiJets_mass.begin(), AK4PuppiJets_mass.end());

    	AK8PuppiJets_pt_float.assign(AK8PuppiJets_pt.begin(), AK8PuppiJets_pt.end()); 
	    AK8PuppiJets_eta_float.assign(AK8PuppiJets_eta.begin(), AK8PuppiJets_eta.end()); 
        AK8PuppiJets_phi_float.assign(AK8PuppiJets_phi.begin(), AK8PuppiJets_phi.end()); 
        AK8PuppiJets_mass_float.assign(AK8PuppiJets_mass.begin(), AK8PuppiJets_mass.end());

        HLTJet80_pt_float.assign(HLTJet80_pt.begin(), HLTJet80_pt.end()); 
        HLTJet80_eta_float.assign(HLTJet80_eta.begin(), HLTJet80_eta.end()); 
        HLTJet80_phi_float.assign(HLTJet80_phi.begin(), HLTJet80_phi.end()); 

        HLTJet60_pt_float.assign(HLTJet60_pt.begin(), HLTJet60_pt.end()); 
        HLTJet60_eta_float.assign(HLTJet60_eta.begin(), HLTJet60_eta.end()); 
        HLTJet60_phi_float.assign(HLTJet60_phi.begin(), HLTJet60_phi.end()); 

        HLTJet80_MatchedCalo_pt_float.assign(HLTJet80_MatchedCalo_pt.begin(), HLTJet80_MatchedCalo_pt.end()); 
        HLTJet80_MatchedCalo_eta_float.assign(HLTJet80_MatchedCalo_eta.begin(), HLTJet80_MatchedCalo_eta.end()); 
        HLTJet80_MatchedCalo_phi_float.assign(HLTJet80_MatchedCalo_phi.begin(), HLTJet80_MatchedCalo_phi.end()); 

        HLTJet60_MatchedCalo_pt_float.assign(HLTJet60_MatchedCalo_pt.begin(), HLTJet60_MatchedCalo_pt.end()); 
        HLTJet60_MatchedCalo_eta_float.assign(HLTJet60_MatchedCalo_eta.begin(), HLTJet60_MatchedCalo_eta.end()); 
        HLTJet60_MatchedCalo_phi_float.assign(HLTJet60_MatchedCalo_phi.begin(), HLTJet60_MatchedCalo_phi.end()); 

        HLTAK4PFJetLoose_pt_float.assign(HLTAK4PFJetLoose_pt.begin(), HLTAK4PFJetLoose_pt.end());
        HLTAK4PFJetLoose_eta_float.assign(HLTAK4PFJetLoose_eta.begin(), HLTAK4PFJetLoose_eta.end());
        HLTAK4PFJetLoose_phi_float.assign(HLTAK4PFJetLoose_phi.begin(), HLTAK4PFJetLoose_phi.end());


        HLTAK4PFJetTight_pt_float.assign(HLTAK4PFJetTight_pt.begin(), HLTAK4PFJetTight_pt.end());
        HLTAK4PFJetTight_eta_float.assign(HLTAK4PFJetTight_eta.begin(), HLTAK4PFJetTight_eta.end());
        HLTAK4PFJetTight_phi_float.assign(HLTAK4PFJetTight_phi.begin(), HLTAK4PFJetTight_phi.end());

        HLTAK4PFJet_pt_float.assign(HLTAK4PFJet_pt.begin(), HLTAK4PFJet_pt.end());
        HLTAK4PFJet_eta_float.assign(HLTAK4PFJet_eta.begin(), HLTAK4PFJet_eta.end());
        HLTAK4PFJet_phi_float.assign(HLTAK4PFJet_phi.begin(), HLTAK4PFJet_phi.end());
//   if(iEvent.id().event() > 709310)
// 	std::cout<<"PIPPO\t before filling 11\n";				                                  
              
        //if (!isMC && passedFullSelection==true) passedEventsTree_All->Fill();        
        if (!isMC && passedTrig==true) passedEventsTree_All->Fill();        
    }    //primary vertex,notDuplicate
    else { if (verbose) cout<<Run<<":"<<LumiSect<<":"<<Event<<" failed primary vertex"<<endl;}
    
    GENjet_pt_float.clear(); GENjet_pt_float.assign(GENjet_pt.begin(),GENjet_pt.end());
    GENjet_eta_float.clear(); GENjet_eta_float.assign(GENjet_eta.begin(),GENjet_eta.end());
    GENjet_phi_float.clear(); GENjet_phi_float.assign(GENjet_phi.begin(),GENjet_phi.end());
    GENjet_mass_float.clear(); GENjet_mass_float.assign(GENjet_mass.begin(),GENjet_mass.end());
 
    //if (isMC && passedFullSelection==true) passedEventsTree_All->Fill();
    if (isMC ) passedEventsTree_All->Fill();
    
    if (nEventsTotal==1000.0) passedEventsTree_All->OptimizeBaskets();
    
}



// ------------ method called once each job just before starting event loop  ------------
void 
HccAna::beginJob()
{
    using namespace edm;
    using namespace std;
    using namespace pat;
		
		
    bookPassedEventTree("passedEvents", passedEventsTree_All);
		
    firstEntry = true;

}

// ------------ method called once each job just after ending the event loop  ------------
void 
HccAna::endJob() 
{
    histContainer_["NEVENTS"]->SetBinContent(1,nEventsTotal);
    histContainer_["NEVENTS"]->GetXaxis()->SetBinLabel(1,"N Events in Sample");
    histContainer_["SUMWEIGHTS"]->SetBinContent(1,sumWeightsTotal);
    histContainer_["SUMWEIGHTSPU"]->SetBinContent(1,sumWeightsTotalPU);
    histContainer_["SUMWEIGHTS"]->GetXaxis()->SetBinLabel(1,"sum Weights in Sample");
    histContainer_["SUMWEIGHTSPU"]->GetXaxis()->SetBinLabel(1,"sum Weights PU in Sample");
}

void
HccAna::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{

    //massErr.init(iSetup);
    if (isMC) {
        edm::Handle<LHERunInfoProduct> run;
        typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
        try {

            int pos=0;
            iRun.getByLabel( edm::InputTag("externalLHEProducer"), run );
            LHERunInfoProduct myLHERunInfoProduct = *(run.product());
            typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
            for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
                std::cout << iter->tag() << std::endl;
                std::vector<std::string> lines = iter->lines();
                for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
                    std::string pdfid=lines.at(iLine);
                    if (pdfid.substr(1,6)=="weight" && pdfid.substr(8,2)=="id") {
                        std::cout<<pdfid<<std::endl;
                        std::string pdf_weight_id = pdfid.substr(12,4);
                        int pdf_weightid=atoi(pdf_weight_id.c_str());
//                         std::cout<<"parsed id: "<<pdf_weightid<<std::endl;
                        if (pdf_weightid==2001) {posNNPDF=int(pos);}
                        pos+=1;
                    }
                }
            }
        }
        catch(...) {
            std::cout<<"No LHERunInfoProduct"<<std::endl;
        }
    }

}


// ------------ method called when ending the processing of a run  ------------
void 
HccAna::endRun(const edm::Run& iRun, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
HccAna::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
HccAna::endLuminosityBlock(edm::LuminosityBlock const& lumiSeg,edm::EventSetup const& eSetup)
{
    using namespace edm;
    using namespace std;
    // Keep track of all the events run over
    edm::Handle<MergeableCounter> numEventsCounter;
    lumiSeg.getByLabel("nEventsTotal", numEventsCounter);    
    if(numEventsCounter.isValid()) {
        std::cout<<"numEventsCounter->value "<<numEventsCounter->value<<endl;
        nEventsTotal += numEventsCounter->value;        
    }
}

// ============================ UF Functions ============================= //



void HccAna::bookPassedEventTree(TString treeName, TTree *tree)
{     


    using namespace edm;
    using namespace pat;
    using namespace std;

    // -------------------------                                                                                                                                                                        
    // RECO level information                                                                                                                                                                           
    // -------------------------                                                                                                                                                                        
    // Event variables
    tree->Branch("Run",&Run,"Run/l");
    tree->Branch("Event",&Event,"Event/l");
    tree->Branch("LumiSect",&LumiSect,"LumiSect/l");
    tree->Branch("nVtx",&nVtx,"nVtx/I");
    tree->Branch("nInt",&nInt,"nInt/I");
    tree->Branch("puN", &puN, "puN/I");
    tree->Branch("PV_x", &PV_x, "PV_x/F");
    tree->Branch("PV_y", &PV_y, "PV_y/F");
    tree->Branch("PV_z", &PV_z, "PV_z/F");
    tree->Branch("BS_x", &BS_x, "BS_x/F");
    tree->Branch("BS_y", &BS_y, "BS_y/F");
    tree->Branch("BS_z", &BS_z, "BS_z/F");
    tree->Branch("BS_xErr", &BS_xErr, "BS_xErr/F");
    tree->Branch("BS_yErr", &BS_yErr, "BS_yErr/F");
    tree->Branch("BS_zErr", &BS_zErr, "BS_zErr/F");
    tree->Branch("BeamWidth_x", &BeamWidth_x, "BeamWidth_x/F");
    tree->Branch("BeamWidth_y", &BeamWidth_y, "BeamWidth_y/F");
    tree->Branch("BeamWidth_xErr", &BeamWidth_xErr, "BeamWidth_xErr/F");
    tree->Branch("BeamWidth_yErr", &BeamWidth_yErr, "BeamWidth_yErr/F");
    tree->Branch("finalState",&finalState,"finalState/I");
    tree->Branch("genWeight",&genWeight);
    tree->Branch("triggersPassed",&triggersPassed);
    tree->Branch("passedTrig",&passedTrig,"passedTrig/O");
    tree->Branch("passedZqqSelection",&passedZqqSelection);
    tree->Branch("Trigger_l1name",&Trigger_l1name);
    tree->Branch("Trigger_l1decision",&Trigger_l1decision);
    tree->Branch("Trigger_hltname",&Trigger_hltname);
    tree->Branch("Trigger_hltdecision",&Trigger_hltdecision);
		

    //tree->Branch("ALLlep_id",&ALLlep_id);
    //tree->Branch("ALLlep_pt",&ALLlep_pt);
    //tree->Branch("ALLlep_eta",&ALLlep_eta);
    //tree->Branch("ALLlep_phi",&ALLlep_phi);
    //tree->Branch("ALLlep_mass",&ALLlep_mass);
    tree->Branch("Ele_id",&Ele_id);
    tree->Branch("Ele_pt",&Ele_pt);
    tree->Branch("Ele_isPassID",&Ele_isPassID);
    tree->Branch("Ele_eta",&Ele_eta);
    tree->Branch("Ele_phi",&Ele_phi);
    tree->Branch("Ele_mass",&Ele_mass);
    tree->Branch("Ele_dxy",&Ele_dxy);
    tree->Branch("Ele_dz",&Ele_dz);
    tree->Branch("Ele_hcalIso", &Ele_hcalIso);
    tree->Branch("Ele_ecalIso", &Ele_ecalIso); 
    tree->Branch("Ele_trackIso", &Ele_trackIso);
    tree->Branch("Ele_isEB", &Ele_isEB);
    tree->Branch("Ele_IsoCal", &Ele_IsoCal);
    tree->Branch("Ele_03_Neutral", &Ele_03_Neutral);

    //tree->Branch("Ele_PF_Iso_R04",&Ele_PF_Iso_R04);
    tree->Branch("Muon_id",&Muon_id);
    tree->Branch("Muon_pt",&Muon_pt);
    tree->Branch("Muon_eta",&Muon_eta);
    tree->Branch("Muon_phi",&Muon_phi);
    tree->Branch("Muon_PF_Iso_R04",&Muon_PF_Iso_R04);
    tree->Branch("Muon_mass",&Muon_mass);
    tree->Branch("Muon_dxy",&Muon_dxy);
    tree->Branch("Muon_dz",&Muon_dz);
    tree->Branch("Muon_PassLooseID", &Muon_PassLooseID);
    tree->Branch("Muon_isPF",&Muon_isPF);
    tree->Branch("AK4lep_id",&AK4lep_id);
    tree->Branch("AK4lep_pt",&AK4lep_pt);
    tree->Branch("AK4lep_eta",&AK4lep_eta);
    tree->Branch("AK4lep_phi",&AK4lep_phi);
    tree->Branch("AK4lep_mass",&AK4lep_mass);


    // MET
    tree->Branch("met",&met,"met/F");
    tree->Branch("met_pt",&met_pt,"met_pt/F");
    tree->Branch("met_phi",&met_phi,"met_phi/F");
    /*tree->Branch("met_jesup",&met_jesup,"met_jesup/F");
    tree->Branch("met_phi_jesup",&met_phi_jesup,"met_phi_jesup/F");
    tree->Branch("met_jesdn",&met_jesdn,"met_jesdn/F");
    tree->Branch("met_phi_jesdn",&met_phi_jesdn,"met_phi_jesdn/F");
    tree->Branch("met_uncenup",&met_uncenup,"met_uncenup/F");
    tree->Branch("met_phi_uncenup",&met_phi_uncenup,"met_phi_uncenup/F");
    tree->Branch("met_uncendn",&met_uncendn,"met_uncendn/F");
    tree->Branch("met_phi_uncendn",&met_phi_uncendn,"met_phi_uncendn/F");*/

    // Jets
    tree->Branch("n_jets", &n_jets);

    // Puppi AK4jets with ParticleNet taggers
    tree->Branch("AK4PuppiJets_pt",&AK4PuppiJets_pt_float);
    tree->Branch("AK4PuppiJetsJESUp_pt",&AK4PuppiJetsJESUp_pt_float);
    tree->Branch("AK4PuppiJetsJESDown_pt",&AK4PuppiJetsJESDown_pt_float);
    tree->Branch("AK4PuppiJetsSmear_pt",&AK4PuppiJetsSmear_pt_float);
    tree->Branch("AK4PuppiJetsSmearUp_pt",&AK4PuppiJetsSmearUp_pt_float);
    tree->Branch("AK4PuppiJetsSmearDown_pt",&AK4PuppiJetsSmearDown_pt_float);
    tree->Branch("AK4PuppiJets_uncorrpt",&AK4PuppiJets_uncorrpt_float);
    tree->Branch("AK4PuppiJets_eta",&AK4PuppiJets_eta_float);
    tree->Branch("AK4PuppiJets_phi",&AK4PuppiJets_phi_float);
    tree->Branch("AK4PuppiJets_mass",&AK4PuppiJets_mass_float);   
    tree->Branch("AK4PuppiJets_hadrFlav",&AK4PuppiJets_hadrFlav);   
    tree->Branch("AK4PuppiJets_partFlav",&AK4PuppiJets_partFlav);   
    tree->Branch("AK4PuppiJets_qgl",&AK4PuppiJets_qgl);   


    //ParticleNet discriminants
    tree->Branch("jet_pfParticleNetAK4JetTags_probb", &jet_pfParticleNetAK4JetTags_probb);	
    tree->Branch("jet_pfParticleNetAK4JetTags_probc", &jet_pfParticleNetAK4JetTags_probc);	
    tree->Branch("jet_pfParticleNetAK4JetTags_probuds", &jet_pfParticleNetAK4JetTags_probuds);	
    tree->Branch("jet_pfParticleNetAK4JetTags_probg", &jet_pfParticleNetAK4JetTags_probg);	
    tree->Branch("jet_pfParticleNetAK4JetTags_probtauh", &jet_pfParticleNetAK4JetTags_probtauh);
    tree->Branch("jet_PNetRegPtRawCorr", &jet_PNetRegPtRawCorr);	
    tree->Branch("jet_PNetRegPtRawCorrNeutrino", &jet_PNetRegPtRawCorrNeutrino);	
    tree->Branch("jet_PNetRegPtRawRes ", &jet_PNetRegPtRawRes);	

    tree->Branch("jet_pfParticleNetAK4JetTags_CvsB", &jet_pfParticleNetAK4JetTags_CvsB);	
    tree->Branch("jet_pfParticleNetAK4JetTags_CvsL", &jet_pfParticleNetAK4JetTags_CvsL);	
    tree->Branch("jet_pfParticleNetAK4JetTags_CvsAll", &jet_pfParticleNetAK4JetTags_CvsAll);	
    tree->Branch("jet_pfParticleNetAK4JetTags_BvsC", &jet_pfParticleNetAK4JetTags_BvsC);	
    tree->Branch("jet_pfParticleNetAK4JetTags_BvsL", &jet_pfParticleNetAK4JetTags_BvsL);	
    tree->Branch("jet_pfParticleNetAK4JetTags_BvsAll", &jet_pfParticleNetAK4JetTags_BvsAll);
    tree->Branch("jet_pfParticleNetAK4JetTags_QvsG", &jet_pfParticleNetAK4JetTags_QvsG);	
    //DeepJet discriminants
    tree->Branch("jet_pfDeepJetAK4JetTags_probb",&jet_pfDeepJetAK4JetTags_probb);
    tree->Branch("jet_pfDeepJetAK4JetTags_probbb",&jet_pfDeepJetAK4JetTags_probbb);
    tree->Branch("jet_pfDeepJetAK4JetTags_problepb",&jet_pfDeepJetAK4JetTags_problepb); 
    tree->Branch("jet_pfDeepJetAK4JetTags_probc",&jet_pfDeepJetAK4JetTags_probc);
    tree->Branch("jet_pfDeepJetAK4JetTags_probuds",&jet_pfDeepJetAK4JetTags_probuds);
    tree->Branch("jet_pfDeepJetAK4JetTags_probg",&jet_pfDeepJetAK4JetTags_probg);

    //DeepCSV discriminants
    tree->Branch("jet_pfDeepCSVAK4JetTags_probb",&jet_pfDeepCSVAK4JetTags_probb);
    tree->Branch("jet_pfDeepCSVAK4JetTags_probbb",&jet_pfDeepCSVAK4JetTags_probbb);
    tree->Branch("jet_pfDeepCSVAK4JetTags_probc",&jet_pfDeepCSVAK4JetTags_probc);
    tree->Branch("jet_pfDeepCSVAK4JetTags_probudsg",&jet_pfDeepCSVAK4JetTags_probudsg);

	// Puppi AK8jets with ParticleNet(-MD) and DeepDoubleX taggers
    tree->Branch("leadingAK8_pt_idx",&leadingAK8_pt_idx);
    tree->Branch("subleadingAK8_pt_idx",&subleadingAK8_pt_idx);

	tree->Branch("AK8PuppiJets_pt",&AK8PuppiJets_pt_float);	
	tree->Branch("AK8PuppiJets_eta",&AK8PuppiJets_eta_float);
	tree->Branch("AK8PuppiJets_phi",&AK8PuppiJets_phi_float);
	tree->Branch("AK8PuppiJets_mass",&AK8PuppiJets_mass_float);
    tree->Branch("AK8PuppiJets_softdropmass",&AK8PuppiJets_softdropmass);
	
	tree->Branch("jet_pfParticleNetJetTags_probZbb", &jet_pfParticleNetJetTags_probZbb);
	tree->Branch("jet_pfParticleNetJetTags_probZcc", &jet_pfParticleNetJetTags_probZcc);
	tree->Branch("jet_pfParticleNetJetTags_probZqq", &jet_pfParticleNetJetTags_probZqq);
	tree->Branch("jet_pfParticleNetJetTags_probQCDbb", &jet_pfParticleNetJetTags_probQCDbb);
	tree->Branch("jet_pfParticleNetJetTags_probQCDcc", &jet_pfParticleNetJetTags_probQCDcc);
	tree->Branch("jet_pfParticleNetJetTags_probQCDb", &jet_pfParticleNetJetTags_probQCDb);
	tree->Branch("jet_pfParticleNetJetTags_probQCDc", &jet_pfParticleNetJetTags_probQCDc);
	tree->Branch("jet_pfParticleNetJetTags_probQCDothers", &jet_pfParticleNetJetTags_probQCDothers);
	tree->Branch("jet_pfParticleNetJetTags_probHbb", &jet_pfParticleNetJetTags_probHbb);
	tree->Branch("jet_pfParticleNetJetTags_probHcc", &jet_pfParticleNetJetTags_probHcc);
	tree->Branch("jet_pfParticleNetJetTags_probHqqqq", &jet_pfParticleNetJetTags_probHqqqq);

	tree->Branch("jet_pfMassDecorrelatedParticleNetJetTags_probXbb", &jet_pfMassDecorrelatedParticleNetJetTags_probXbb);
    tree->Branch("jet_pfMassDecorrelatedParticleNetJetTags_probXcc", &jet_pfMassDecorrelatedParticleNetJetTags_probXcc);
    tree->Branch("jet_pfMassDecorrelatedParticleNetJetTags_probXqq", &jet_pfMassDecorrelatedParticleNetJetTags_probXqq);
    tree->Branch("jet_pfMassDecorrelatedParticleNetJetTags_probQCDbb", &jet_pfMassDecorrelatedParticleNetJetTags_probQCDbb);
    tree->Branch("jet_pfMassDecorrelatedParticleNetJetTags_probQCDcc", &jet_pfMassDecorrelatedParticleNetJetTags_probQCDcc);
    tree->Branch("jet_pfMassDecorrelatedParticleNetJetTags_probQCDb", &jet_pfMassDecorrelatedParticleNetJetTags_probQCDb);
    tree->Branch("jet_pfMassDecorrelatedParticleNetJetTags_probQCDc", &jet_pfMassDecorrelatedParticleNetJetTags_probQCDc);
    tree->Branch("jet_pfMassDecorrelatedParticleNetJetTags_probQCDothers", &jet_pfMassDecorrelatedParticleNetJetTags_probQCDothers);

	tree->Branch("jet_pfMassIndependentDeepDoubleBvLV2JetTags_probHbb", &jet_pfMassIndependentDeepDoubleBvLV2JetTags_probHbb);
	tree->Branch("jet_pfMassIndependentDeepDoubleCvLV2JetTags_probHcc", &jet_pfMassIndependentDeepDoubleCvLV2JetTags_probHcc);
	tree->Branch("jet_pfMassIndependentDeepDoubleCvBV2JetTags_probHcc", &jet_pfMassIndependentDeepDoubleCvBV2JetTags_probHcc);

    //HLT jets for turn on curves and scale factors
    tree->Branch("HLTJet80_pt",&HLTJet80_pt_float);
    tree->Branch("HLTJet80_eta",&HLTJet80_eta_float);
    tree->Branch("HLTJet80_phi",&HLTJet80_phi_float);

    tree->Branch("HLTJet60_pt",&HLTJet60_pt_float);
    tree->Branch("HLTJet60_eta",&HLTJet60_eta_float);
    tree->Branch("HLTJet60_phi",&HLTJet60_phi_float);

    tree->Branch("HLTJet80_MatchedCalo_pt",&HLTJet80_MatchedCalo_pt_float);
    tree->Branch("HLTJet80_MatchedCalo_eta",&HLTJet80_MatchedCalo_eta_float);
    tree->Branch("HLTJet80_MatchedCalo_phi",&HLTJet80_MatchedCalo_phi_float);

    tree->Branch("HLTJet60_MatchedCalo_pt",&HLTJet60_MatchedCalo_pt_float);
    tree->Branch("HLTJet60_MatchedCalo_eta",&HLTJet60_MatchedCalo_eta_float);
    tree->Branch("HLTJet60_MatchedCalo_phi",&HLTJet60_MatchedCalo_phi_float);

    tree->Branch("HLTAK4PFJetLoose_pt",&HLTAK4PFJetLoose_pt_float);
    tree->Branch("HLTAK4PFJetLoose_eta",&HLTAK4PFJetLoose_eta_float);
    tree->Branch("HLTAK4PFJetLoose_phi",&HLTAK4PFJetLoose_phi_float);

    tree->Branch("HLTAK4PFJetTight_pt",&HLTAK4PFJetTight_pt_float);
    tree->Branch("HLTAK4PFJetTight_eta",&HLTAK4PFJetTight_eta_float);
    tree->Branch("HLTAK4PFJetTight_phi",&HLTAK4PFJetTight_phi_float);

    tree->Branch("HLTAK4PFJet_pt",&HLTAK4PFJet_pt_float);
    tree->Branch("HLTAK4PFJet_eta",&HLTAK4PFJet_eta_float);
    tree->Branch("HLTAK4PFJet_phi",&HLTAK4PFJet_phi_float);

    //hlt jets
    tree->Branch("hltjetForBTag_pt",&hltjetForBTag_pt_float);
    tree->Branch("hltjetForBTag_eta",&hltjetForBTag_eta_float);
    tree->Branch("hltjetForBTag_phi",&hltjetForBTag_phi_float);
    tree->Branch("hltjetForBTag_mass",&hltjetForBTag_mass_float);
    tree->Branch("hltjetForBTag_ParticleNet_probb",&hltParticleNetONNXJetTags_probb);
    tree->Branch("hltjetForBTag_ParticleNet_probc",&hltParticleNetONNXJetTags_probc);
    tree->Branch("hltjetForBTag_ParticleNet_probuds",&hltParticleNetONNXJetTags_probuds);
    tree->Branch("hltjetForBTag_ParticleNet_probg",&hltParticleNetONNXJetTags_probg);
    tree->Branch("hltjetForBTag_ParticleNet_probtauh",&hltParticleNetONNXJetTags_probtauh);

    //hltAK4PFJetsCorrected
    tree->Branch("hltAK4PFJetsCorrected_pt",&hltAK4PFJetsCorrected_pt_float);
    tree->Branch("hltAK4PFJetsCorrected_eta",&hltAK4PFJetsCorrected_eta_float);
    tree->Branch("hltAK4PFJetsCorrected_phi",&hltAK4PFJetsCorrected_phi_float);
    tree->Branch("hltAK4PFJetsCorrected_mass",&hltAK4PFJetsCorrected_mass_float);

    //L1 jets
    tree->Branch("L1jet_pt",&L1jet_pt_float);
    tree->Branch("L1jet_eta",&L1jet_eta_float);
    tree->Branch("L1jet_phi",&L1jet_phi_float);
    tree->Branch("L1jet_mass",&L1jet_mass_float);
	
    //L1 muons
    tree->Branch("L1muon_pt",&L1muon_pt_float);
    tree->Branch("L1muon_eta",&L1muon_eta_float);
    tree->Branch("L1muon_phi",&L1muon_phi_float);
    tree->Branch("L1muon_mass",&L1muon_mass_float);
    tree->Branch("L1muon_qual",&L1muon_qual);
		
    //L1 HT
    tree->Branch("L1ht",&L1ht, "L1ht/F");
	

    // Event Category
    tree->Branch("EventCat",&EventCat,"EventCat/I");

    // -------------------------                                                                                                                                                                        
    // GEN level information                                                                                                                                                                            
    // -------------------------                                                                                                                                                                        
    //Event variables
    tree->Branch("GENfinalState",&GENfinalState,"GENfinalState/I");

    //quark
    tree->Branch("quark_pt", &quark_pt_float);
    tree->Branch("quark_eta", &quark_eta_float);
    tree->Branch("quark_phi", &quark_phi_float);
    tree->Branch("quark_flavour", &quark_flavour);
    tree->Branch("quark_VBF", &quark_VBF);
    
    // Z boson
    tree->Branch("Z_pt", &Z_pt);
    tree->Branch("Z_eta", &Z_eta);
    tree->Branch("Z_phi", &Z_phi);
    tree->Branch("Z_mass", &Z_mass);

    // Jets
    tree->Branch("n_GENjets", &n_GENjets);
    tree->Branch("GENjet_pt",&GENjet_pt_float);
    tree->Branch("GENjet_eta",&GENjet_eta_float);
    tree->Branch("GENjet_phi",&GENjet_phi_float);
    tree->Branch("GENjet_mass",&GENjet_mass_float);

}

/*void HccAna::setTreeVariables( const edm::Event& iEvent, const edm::EventSetup& iSetup,
                                   std::vector<pat::Muon> selectedMuons, std::vector<pat::Electron> selectedElectrons, 
                                   std::vector<pat::Muon> recoMuons, std::vector<pat::Electron> recoElectrons, 
                                   std::vector<pat::Jet> goodJets, std::vector<float> goodJetQGTagger, 
                                   std::vector<float> goodJetaxis2, std::vector<float> goodJetptD, std::vector<int> goodJetmult,
                                   std::vector<pat::Jet> selectedMergedJets,
                                   std::map<unsigned int, TLorentzVector> selectedFsrMap)*/
void HccAna::setTreeVariables( const edm::Event& iEvent, const edm::EventSetup& iSetup,
                                   //std::vector<pat::Jet> goodJets, //std::vector<float> goodJetQGTagger,
                                   //std::vector<float> goodJetaxis2, std::vector<float> goodJetptD, std::vector<int> goodJetmult,
                                   //std::vector<pat::Jet> selectedMergedJets,
                                   edm::Handle<edm::View<pat::Jet> > AK4PuppiJets,
                                   //edm::Handle<edm::View<pat::Jet> > AK4PuppiJetsJESUp,
                                   //edm::Handle<edm::View<pat::Jet> > AK4PuppiJetsJESDown,
                                   edm::Handle<edm::View<pat::Jet> > AK4PuppiJetsUncorr,
                                   JetCorrectionUncertainty *uncAK4,
                                   edm::Handle<edm::View<pat::Jet> > AK4PuppiJetsSmear,
                                   edm::Handle<edm::View<pat::Jet> > AK4PuppiJetsSmearUp,
                                   edm::Handle<edm::View<pat::Jet> > AK4PuppiJetsSmearDown,
                                   edm::Handle<edm::View<pat::Jet> > AK8PuppiJets,
                                  // const edm::TriggerNames trigNames,
                                 // edm::Handle<pat::TriggerObjectStandAlone> triggerObjects,
                                 //edm::Handle<std::vector<reco::PFJet>> hltjets,
                                 //edm::Handle<edm::View<reco::PFJet>> hltjetsForBTag,
                                 //edm::Handle<edm::View<reco::PFJet>> hltAK4PFJetsCorrected,
                                 //edm::Handle<reco::JetTagCollection> pfJetTagCollectionParticleNetprobc,
                                 //edm::Handle<reco::JetTagCollection> pfJetTagCollectionParticleNetprobb,
                                 //edm::Handle<reco::JetTagCollection> pfJetTagCollectionParticleNetprobuds,
                                 //edm::Handle<reco::JetTagCollection> pfJetTagCollectionParticleNetprobg,
                                 //edm::Handle<reco::JetTagCollection> pfJetTagCollectionParticleNetprobtauh,
                                 //  edm::Handle<BXVector<l1t::Jet> > bxvCaloJets,
                                 //  edm::Handle<BXVector<l1t::Muon> > bxvCaloMuons,
                                  // edm::Handle<BXVector<l1t::EtSum> > bxvCaloHT,
                                 //edm::Handle<edm::View<pat::Muon> > muons,
                                 //edm::Handle<edm::View<pat::Electron> > electrons)
                                   std::vector<pat::Muon> AllMuons, std::vector<pat::Electron> AllElectrons, const reco::Vertex *ver)
{

   

    using namespace edm;
    using namespace pat;
    using namespace std;

	for(unsigned int jmu=0; jmu<AllMuons.size(); jmu++){
        /*ALLlep_pt.push_back(AllMuons[jmu].pt());
        ALLlep_eta.push_back(AllMuons[jmu].eta());
        ALLlep_phi.push_back(AllMuons[jmu].phi());
        ALLlep_mass.push_back(AllMuons[jmu].mass());
        ALLlep_id.push_back(AllMuons[jmu].pdgId());*/
        Muon_pt.push_back(AllMuons[jmu].pt());
        Muon_eta.push_back(AllMuons[jmu].eta());
        Muon_phi.push_back(AllMuons[jmu].phi());
        Muon_mass.push_back(AllMuons[jmu].mass());
        Muon_dxy.push_back(AllMuons[jmu].innerTrack()->dxy(ver->position()));
        Muon_dz.push_back(AllMuons[jmu].innerTrack()->dz(ver->position()));
        Muon_id.push_back(AllMuons[jmu].pdgId());
        Muon_isPF.push_back(AllMuons[jmu].isPFMuon());
        Muon_PF_Iso_R04.push_back((AllMuons[jmu].pfIsolationR04().sumChargedHadronPt + TMath::Max(AllMuons[jmu].pfIsolationR04().sumNeutralHadronEt + AllMuons[jmu].pfIsolationR04().sumPhotonEt - AllMuons[jmu].pfIsolationR04().sumPUPt/2.0,0.0))/AllMuons[jmu].pt());
        Muon_PassLooseID.push_back(AllMuons[jmu].isLooseMuon());
//      if(AllMuons[jmu].pt()>20 && abs(AllMuons[jmu].eta())<2.4 && AllMuons[jmu].isLooseMuon() && Mu_PF_Iso_R04<0.4 ){Nmu=Nmu+1;}
        }

	for(unsigned int jel=0; jel<AllElectrons.size(); jel++){
        /*ALLlep_pt.push_back(AllElectrons[jel].pt());
        ALLlep_eta.push_back(AllElectrons[jel].eta());
        ALLlep_phi.push_back(AllElectrons[jel].phi());
        ALLlep_mass.push_back(AllElectrons[jel].mass());
        ALLlep_id.push_back(AllElectrons[jel].pdgId());*/
        Ele_pt.push_back(AllElectrons[jel].pt());
        Ele_eta.push_back(AllElectrons[jel].eta());
        Ele_phi.push_back(AllElectrons[jel].phi());
        Ele_mass.push_back(AllElectrons[jel].mass());
        Ele_dxy.push_back(AllElectrons[jel].gsfTrack()->dxy(ver->position()));
        Ele_dz.push_back(AllElectrons[jel].gsfTrack()->dz(ver->position()));
        Ele_hcalIso.push_back(AllElectrons[jel].dr03TkSumPt() / AllElectrons[jel].pt());
        Ele_ecalIso.push_back(AllElectrons[jel].dr03EcalRecHitSumEt() / AllElectrons[jel].pt()); 
        Ele_trackIso.push_back(AllElectrons[jel].dr03HcalTowerSumEt()/ AllElectrons[jel].pt()); 
        Ele_isEB.push_back(AllElectrons[jel].isEB()); 
        if(AllElectrons[jel].isEB()==true){
            Ele_IsoCal.push_back((std::max(0., AllElectrons[jel].dr03EcalRecHitSumEt() - 1.) + AllElectrons[jel].dr03HcalTowerSumEt()) /AllElectrons[jel].pt());
        } //-1 -> pedestal subtraction needed only in the barrel
        else{
            Ele_IsoCal.push_back((AllElectrons[jel].dr03EcalRecHitSumEt() + AllElectrons[jel].dr03HcalTowerSumEt()) /AllElectrons[jel].pt());
        }
        double El_03_Neutral = AllElectrons[jel].ecalPFClusterIso();
        El_03_Neutral += AllElectrons[jel].hcalPFClusterIso();
        Ele_03_Neutral.push_back(El_03_Neutral);

        Ele_id.push_back(AllElectrons[jel].pdgId());
         //Ele_PF_Iso_R04.push_back((AllElectrons[jel].setIsolation04().sumChargedHadronPt + TMath::Max(AllElectrons[jel].setIsolation04().sumNeutralHadronEt + AllElectrons[jel].setIsolation04().sumPhotonEt - AllElectrons[jel].setIsolation04().sumPUPt/2.0,0.0))/AllElectrons[jel].pt());
        Ele_isPassID.push_back(AllElectrons[jel].electronID("cutBasedElectronID-Fall17-94X-V2-veto"));
//                if(AllElectrons[jel].pt()>20 && abs(AllElectrons[jel].eta())<2.4 && isPassID){Ne=Ne+1;}
        }
       	
    //L1 jets Variables
    /*for (std::vector<l1t::Jet>::const_iterator l1jet = bxvCaloJets->begin(0); l1jet != bxvCaloJets->end(0); ++l1jet) {
      L1jet_pt.push_back(l1jet->pt());
      L1jet_eta.push_back(l1jet->eta());
      L1jet_phi.push_back(l1jet->phi());
      L1jet_mass.push_back(l1jet->mass());
    }

    //L1 muon Variables
    for (std::vector<l1t::Muon>::const_iterator l1muon = bxvCaloMuons->begin(0); l1muon != bxvCaloMuons->end(0); ++l1muon) {
      L1muon_pt.push_back(l1muon->pt());
      L1muon_eta.push_back(l1muon->eta());
      L1muon_phi.push_back(l1muon->phi());
      L1muon_mass.push_back(l1muon->mass());
      L1muon_qual.push_back(l1muon->hwQual());
    }

    //L1 HT sum
    for (std::vector<l1t::EtSum>::const_iterator l1Et = bxvCaloHT->begin(0); l1Et != bxvCaloHT->end(0); ++l1Et) {
      if (l1Et->getType() == l1t::EtSum::EtSumType::kTotalHt){
        L1ht= l1Et->et();
      }
    }*/


	
    //hltAK4PFJetsCorrected
    /*for(unsigned int ijet=0; ijet<hltAK4PFJetsCorrected->size(); ijet++){
	  //std::cout<<"index jet: "<<ijet<<std::endl;
      //std::cout<<"jet pt: "<<hltjets->at(ijet).pt()<<std::endl;
	  hltAK4PFJetsCorrected_pt.push_back(hltAK4PFJetsCorrected->at(ijet).pt());
	  hltAK4PFJetsCorrected_eta.push_back(hltAK4PFJetsCorrected->at(ijet).eta());
	  hltAK4PFJetsCorrected_phi.push_back(hltAK4PFJetsCorrected->at(ijet).phi());
	  hltAK4PFJetsCorrected_mass.push_back(hltAK4PFJetsCorrected->at(ijet).mass());
    }*/    
    
    //Puppi AK8jets with ParticleNet and DeepDoubleX taggers    
    
    float subleadingAK8_pt=-1000;
    float leadingAK8_pt= -100;
 
    for(unsigned int jjet=0; jjet<AK8PuppiJets->size(); jjet++){
      AK8PuppiJets_pt.push_back(AK8PuppiJets->at(jjet).pt());
      AK8PuppiJets_eta.push_back(AK8PuppiJets->at(jjet).eta());
      AK8PuppiJets_phi.push_back(AK8PuppiJets->at(jjet).phi());
      AK8PuppiJets_mass.push_back(AK8PuppiJets->at(jjet).mass());
      AK8PuppiJets_softdropmass.push_back(AK8PuppiJets->at(jjet).userFloat("ak8PFJetsPuppiSoftDropMass"));

      if(AK8PuppiJets->at(jjet).pt()>leadingAK8_pt){
        subleadingAK8_pt = leadingAK8_pt;
        leadingAK8_pt=AK8PuppiJets->at(jjet).pt();
        subleadingAK8_pt_idx = leadingAK8_pt_idx;
        leadingAK8_pt_idx=jjet;}

      if(AK8PuppiJets->at(jjet).pt()<leadingAK8_pt && AK8PuppiJets->at(jjet).pt()>subleadingAK8_pt){
        subleadingAK8_pt = AK8PuppiJets->at(jjet).pt();
        subleadingAK8_pt_idx = jjet;}

      jet_pfParticleNetJetTags_probZbb.push_back(AK8PuppiJets->at(jjet).bDiscriminator("pfParticleNetJetTags:probZbb"));
      jet_pfParticleNetJetTags_probZcc.push_back(AK8PuppiJets->at(jjet).bDiscriminator("pfParticleNetJetTags:probZcc"));
      jet_pfParticleNetJetTags_probZqq.push_back(AK8PuppiJets->at(jjet).bDiscriminator("pfParticleNetJetTags:probZqq"));
      jet_pfParticleNetJetTags_probQCDbb.push_back(AK8PuppiJets->at(jjet).bDiscriminator("pfParticleNetJetTags:probQCDbb"));
      jet_pfParticleNetJetTags_probQCDcc.push_back(AK8PuppiJets->at(jjet).bDiscriminator("pfParticleNetJetTags:probQCDcc"));
      jet_pfParticleNetJetTags_probQCDb.push_back(AK8PuppiJets->at(jjet).bDiscriminator("pfParticleNetJetTags:probQCDb"));
      jet_pfParticleNetJetTags_probQCDc.push_back(AK8PuppiJets->at(jjet).bDiscriminator("pfParticleNetJetTags:probQCDc"));
      jet_pfParticleNetJetTags_probQCDothers.push_back(AK8PuppiJets->at(jjet).bDiscriminator("pfParticleNetJetTags:probQCDothers"));
      jet_pfParticleNetJetTags_probHbb.push_back(AK8PuppiJets->at(jjet).bDiscriminator("pfParticleNetJetTags:probHbb"));
      jet_pfParticleNetJetTags_probHcc.push_back(AK8PuppiJets->at(jjet).bDiscriminator("pfParticleNetJetTags:probHcc"));
      jet_pfParticleNetJetTags_probHqqqq.push_back(AK8PuppiJets->at(jjet).bDiscriminator("pfParticleNetJetTags:probHqqqq"));
      
      jet_pfMassDecorrelatedParticleNetJetTags_probXbb.push_back(AK8PuppiJets->at(jjet).bDiscriminator("pfMassDecorrelatedParticleNetJetTags:probXbb"));
      jet_pfMassDecorrelatedParticleNetJetTags_probXcc.push_back(AK8PuppiJets->at(jjet).bDiscriminator("pfMassDecorrelatedParticleNetJetTags:probXcc"));
      //std::cout<<"Prob Xbb: "<<AK8PuppiJets->at(jjet).bDiscriminator("pfParticleNetFromMiniAODAK8JetTags:probXbb")<<std::endl;
    	jet_pfMassDecorrelatedParticleNetJetTags_probXqq.push_back(AK8PuppiJets->at(jjet).bDiscriminator("pfMassDecorrelatedParticleNetJetTags:probXqq"));
      jet_pfMassDecorrelatedParticleNetJetTags_probQCDbb.push_back(AK8PuppiJets->at(jjet).bDiscriminator("pfMassDecorrelatedParticleNetJetTags:probQCDbb"));
      jet_pfMassDecorrelatedParticleNetJetTags_probQCDcc.push_back(AK8PuppiJets->at(jjet).bDiscriminator("pfMassDecorrelatedParticleNetJetTags:probQCDcc"));
      jet_pfMassDecorrelatedParticleNetJetTags_probQCDb.push_back(AK8PuppiJets->at(jjet).bDiscriminator("pfMassDecorrelatedParticleNetJetTags:probQCDb"));
      jet_pfMassDecorrelatedParticleNetJetTags_probQCDc.push_back(AK8PuppiJets->at(jjet).bDiscriminator("pfMassDecorrelatedParticleNetJetTags:probQCDc"));
      jet_pfMassDecorrelatedParticleNetJetTags_probQCDothers.push_back(AK8PuppiJets->at(jjet).bDiscriminator("pfMassDecorrelatedParticleNetJetTags:probQCDothers"));
      
      jet_pfMassIndependentDeepDoubleBvLV2JetTags_probHbb.push_back(AK8PuppiJets->at(jjet).bDiscriminator("pfMassIndependentDeepDoubleBvLV2JetTags:probHbb"));// DeepDoubleX discriminator (mass-decorrelation) for H(Z)->bb vs QCD
      jet_pfMassIndependentDeepDoubleCvLV2JetTags_probHcc.push_back(AK8PuppiJets->at(jjet).bDiscriminator("pfMassIndependentDeepDoubleCvLV2JetTags:probHcc"));// DeepDoubleX discriminator (mass-decorrelation) for H(Z)->cc vs QCD
      jet_pfMassIndependentDeepDoubleCvBV2JetTags_probHcc.push_back(AK8PuppiJets->at(jjet).bDiscriminator("pfMassIndependentDeepDoubleCvBV2JetTags:probHcc"));// DeepDoubleX discriminator (mass-decorrelation) for H(Z)->cc vs for H(Z)->bb
      
    }

    //Uncorrected AK4jets puppi

    //for(unsigned int ijet=0; ijet<AK4PuppiJetsUncorr->size(); ijet++){

      // Tight JetID function defined at ../interface/HccJets.h
      //if(jetHelper.patjetID(AK4PuppiJetsUncorr->at(ijet), year)!=2) continue;

    //  AK4PuppiJets_uncorrpt.push_back(AK4PuppiJetsUncorr->at(ijet).pt());
   // }
     
    //Puppi AK4jets with ParticleNet taggers

    //JetCorrectorParameters *pAK4 = new JetCorrectorParameters(JECUncFileAK4.c_str(), "Total") ; 
    //JetCorrectionUncertainty *uncAK4 = new JetCorrectionUncertainty(*pAK4);
    for(unsigned int ijet=0; ijet<AK4PuppiJets->size(); ijet++){

      // Tight JetID function defined at ../interface/HccJets.h
      if(jetHelper.patjetID(AK4PuppiJets->at(ijet), year)!=2) continue;

      AK4PuppiJets_pt.push_back(AK4PuppiJets->at(ijet).pt());
      //cout<<" jet pt: "<<AK4PuppiJets->at(ijet).pt()<<"   jet eta: "<<AK4PuppiJets->at(ijet).eta()<<endl;
      AK4PuppiJets_uncorrpt.push_back(AK4PuppiJetsUncorr->at(ijet).pt()); //uncorrected pt 
      // get jet energy scale Up and Down uncertainties 
      //JetCorrectorParameters *pAK4 = new JetCorrectorParameters(JECUncFileAK4.c_str(), "Total") ; 
      //JetCorrectionUncertainty *uncAK4 = new JetCorrectionUncertainty(*pAK4);
      uncAK4->setJetPt(AK4PuppiJets->at(ijet).pt());
      uncAK4->setJetEta(AK4PuppiJets->at(ijet).eta());
      double vr = uncAK4->getUncertainty(true);
      double vr_up = 1 + vr; //  variation
      //cout<<"vr_up: "<<vr_up<<endl;
      double jesUp= AK4PuppiJets->at(ijet).pt() * vr_up;
      //cout<<"JESUp pt: "<<jesUp<<endl;
      double vr_down = 1 - vr; //  variation
      //cout<<"vr_down: "<<vr_down<<endl;
      double jesDown= AK4PuppiJets->at(ijet).pt() * vr_down;
	  //cout<<"JESDown pt: "<<jesDown<<endl;
      AK4PuppiJetsJESUp_pt.push_back(jesUp);
      AK4PuppiJetsJESDown_pt.push_back(jesDown);
      // jet resolution smearing with Up and Down uncertainties
      AK4PuppiJetsSmear_pt.push_back(AK4PuppiJetsSmear->at(ijet).pt());
      AK4PuppiJetsSmearUp_pt.push_back(AK4PuppiJetsSmearUp->at(ijet).pt());
      AK4PuppiJetsSmearDown_pt.push_back(AK4PuppiJetsSmearDown->at(ijet).pt());
      //cout<<"smearUp pt: "<<AK4PuppiJetsSmearUp->at(ijet).pt()<<"   eta: "<<AK4PuppiJetsSmearUp->at(ijet).eta()<<endl;
      AK4PuppiJets_eta.push_back(AK4PuppiJets->at(ijet).eta());
      AK4PuppiJets_phi.push_back(AK4PuppiJets->at(ijet).phi());
      AK4PuppiJets_mass.push_back(AK4PuppiJets->at(ijet).mass());     
      AK4PuppiJets_hadrFlav.push_back(AK4PuppiJets->at(ijet).hadronFlavour());
      AK4PuppiJets_partFlav.push_back(AK4PuppiJets->at(ijet).partonFlavour());
      //cout<<"parton flavour: "<<AK4PuppiJets->at(ijet).partonFlavour()<<endl; 
      // Tight JetID
      /*double nhf   = AK4PuppiJets->at(ijet).neutralHadronEnergyFraction();
      double nef   = AK4PuppiJets->at(ijet).neutralEmEnergyFraction();
      double chf   = AK4PuppiJets->at(ijet).chargedHadronEnergyFraction();;
      int NumConst = AK4PuppiJets->at(ijet).chargedMultiplicity() + AK4PuppiJets->at(ijet).neutralMultiplicity();
      int CharMult = AK4PuppiJets->at(ijet).chargedMultiplicity();
      int NeutMult = AK4PuppiJets->at(ijet).neutralMultiplicity();*/

      //if (!isTight(AK4PuppiJets->at(ijet).eta(), nhf, nef, chf, NumConst, CharMult, NeutMult)) continue;
      
      float QvsG=-1.;

      if(year==2023){

        jet_pfParticleNetAK4JetTags_probb.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probb"));
        //jet_pfParticleNetAK4JetTags_probb.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetAK4JetTags:probb"));
        jet_pfParticleNetAK4JetTags_probc.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probc"));
        jet_pfParticleNetAK4JetTags_probuds.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probuds"));
        jet_pfParticleNetAK4JetTags_probg.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probg"));
        jet_pfParticleNetAK4JetTags_probtauh.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probtauh"));
        float PNetRegPtRawCorr = AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:ptcorr");
        float PNetRegPtRawCorrNeutrino = AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:ptnu");
        float PNetRegPtRawRes = 0.5*(AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:ptreshigh")-AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:ptreslow"));
        if(fabs(AK4PuppiJets->at(ijet).eta())>=2.5){
          PNetRegPtRawCorr = AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetFromMiniAODAK4PuppiForwardJetTags:ptcorr");
          PNetRegPtRawCorrNeutrino = AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetFromMiniAODAK4PuppiForwardJetTags:ptnu");
          PNetRegPtRawRes = 0.5*(AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetFromMiniAODAK4PuppiForwardJetTags:ptreshigh")-AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetFromMiniAODAK4PuppiForwardJetTags:ptreslow"));
        }
        jet_PNetRegPtRawCorr.push_back(PNetRegPtRawCorr);
        jet_PNetRegPtRawCorrNeutrino.push_back(PNetRegPtRawCorrNeutrino);
        jet_PNetRegPtRawRes.push_back(PNetRegPtRawRes); 
 
        jet_pfParticleNetAK4JetTags_CvsB.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralDiscriminatorsJetTags:CvsB"));
        jet_pfParticleNetAK4JetTags_CvsL.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralDiscriminatorsJetTags:CvsL"));
        jet_pfParticleNetAK4JetTags_CvsAll.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralDiscriminatorsJetTags:CvsAll"));
        //std::cout<<"probc: "<<AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probc")<<std::endl;
        jet_pfParticleNetAK4JetTags_BvsC.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralDiscriminatorsJetTags:BvsC"));
        jet_pfParticleNetAK4JetTags_BvsL.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralDiscriminatorsJetTags:BvsL"));
        jet_pfParticleNetAK4JetTags_BvsAll.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralDiscriminatorsJetTags:BvsAll"));
        QvsG=AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralDiscriminatorsJetTags:QvsG");
        if(fabs(AK4PuppiJets->at(ijet).eta())>=2.5){
          QvsG=AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetFromMiniAODAK4PuppiForwardDiscriminatorsJetTags:QvsG");
        }
        jet_pfParticleNetAK4JetTags_QvsG.push_back(QvsG);
      }

      else{
        jet_pfParticleNetAK4JetTags_probb.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetAK4JetTags:probb"));
        jet_pfParticleNetAK4JetTags_probc.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetAK4JetTags:probc"));
        jet_pfParticleNetAK4JetTags_probuds.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetAK4JetTags:probuds"));
        jet_pfParticleNetAK4JetTags_probg.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetAK4JetTags:probg"));
        jet_pfParticleNetAK4JetTags_probtauh.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetAK4JetTags:probtauh"));

        jet_pfParticleNetAK4JetTags_CvsB.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetAK4DiscriminatorsJetTags:CvsB"));
        jet_pfParticleNetAK4JetTags_CvsL.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetAK4DiscriminatorsJetTags:CvsL"));
        jet_pfParticleNetAK4JetTags_CvsAll.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetAK4DiscriminatorsJetTags:CvsAll"));
        //std::cout<<"probc: "<<AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetAK4DiscriminatorsJetTags:probc")<<std::endl;
        jet_pfParticleNetAK4JetTags_BvsC.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetAK4DiscriminatorsJetTags:BvsC"));
        jet_pfParticleNetAK4JetTags_BvsL.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetAK4DiscriminatorsJetTags:BvsL"));
        jet_pfParticleNetAK4JetTags_BvsAll.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetAK4DiscriminatorsJetTags:BvsAll"));
        QvsG=AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetAK4DiscriminatorsJetTags:QvsG");
        //if(fabs(AK4PuppiJets->at(ijet).eta())>=2.5){
        //  QvsG=AK4PuppiJets->at(ijet).bDiscriminator("pfParticleNetFromMiniAODAK4PuppiForwardJetTags:QvsG");
        //}
        jet_pfParticleNetAK4JetTags_QvsG.push_back(QvsG);

      }

      float qgl= AK4PuppiJets->at(ijet).userFloat("QGTagger:qgLikelihood");
      AK4PuppiJets_qgl.push_back(qgl);
      //cout<<"qgl: "<< qgl<<"     pt: "<<AK4PuppiJets->at(ijet).pt()<<"     eta: "<<AK4PuppiJets->at(ijet).eta()<<endl;

      jet_pfDeepJetAK4JetTags_probb.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfDeepFlavourJetTags:probb")); 
      jet_pfDeepJetAK4JetTags_probbb.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfDeepFlavourJetTags:probbb")); 
      jet_pfDeepJetAK4JetTags_problepb.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfDeepFlavourJetTags:problepb")); 
      jet_pfDeepJetAK4JetTags_probc.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfDeepFlavourJetTags:probc")); 
      jet_pfDeepJetAK4JetTags_probuds.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfDeepFlavourJetTags:probuds")); 
      jet_pfDeepJetAK4JetTags_probg.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfDeepFlavourJetTags:probg"));

      jet_pfDeepCSVAK4JetTags_probb.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfDeepCSVJetTags:probb")); 
      jet_pfDeepCSVAK4JetTags_probbb.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfDeepCSVJetTags:probbb")); 
      jet_pfDeepCSVAK4JetTags_probc.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfDeepCSVJetTags:probc")); 
      jet_pfDeepCSVAK4JetTags_probudsg.push_back(AK4PuppiJets->at(ijet).bDiscriminator("pfDeepCSVJetTags:probudsg"));
      
      /*bool passPFtightID_LepVeto = false;

      float NHF = AK4PuppiJets->at(ijet).neutralHadronEnergyFraction();
      float NEF = AK4PuppiJets->at(ijet).neutralEmEnergyFraction();
      float NumConst = AK4PuppiJets->at(ijet).chargedMultiplicity()+ AK4PuppiJets->at(ijet).neutralMultiplicity();
      float MUF = AK4PuppiJets->at(ijet).muonEnergyFraction();
      float CHF = AK4PuppiJets->at(ijet).chargedHadronEnergyFraction();
      float CHM = AK4PuppiJets->at(ijet).chargedMultiplicity();
      float CEF = AK4PuppiJets->at(ijet).chargedEmEnergyFraction();

      if(NHF < 0.9 && NEF < 0.9 && NumConst > 1 && MUF < 0.8 && CHF > 0.01 && CHM > 0 && CEF <0.8){passPFtightID_LepVeto = true;}

      bool passPUtightID = false;

            if(AK4PuppiJets->at(ijet).pt()>20 && AK4PuppiJets->at(ijet).pt()<30 && AK4PuppiJets->at(ijet).userFloat("pileupJetId:fullDiscriminant")>0.69){passPUtightID = true;}
	      if(AK4PuppiJets->at(ijet).pt()>30 && AK4PuppiJets->at(ijet).pt()<50 && AK4PuppiJets->at(ijet).userFloat("pileupJetId:fullDiscriminant")>0.86){passPUtightID = true;}
      if(AK4PuppiJets->at(ijet).pt()>50){passPUtightID = true;}

      bool overlaps_loose_lepton = false;*/
    }

    for( unsigned int kmu = 0; kmu < AllMuons.size(); kmu++) {
      for(unsigned int kk=0; kk < AK4PuppiJets->size(); kk++){
        bool isMuonFound = false;			
        double this_dR_AKmu = deltaR(AK4PuppiJets->at(kk).eta(), AK4PuppiJets->at(kk).phi(), AllMuons[kmu].eta(), AllMuons[kmu].phi());
          if(this_dR_AKmu<0.4 && !isMuonFound){ 	
            AK4lep_pt.push_back(AllMuons[kmu].pt());
            AK4lep_eta.push_back(AllMuons[kmu].eta());
            AK4lep_phi.push_back(AllMuons[kmu].phi());
            AK4lep_mass.push_back(AllMuons[kmu].mass());
            AK4lep_id.push_back(AllMuons[kmu].pdgId());
            isMuonFound = true;//stop jet cicle if the muon is asscociated to one AK4 jet
          }	
        }
      }
      
      for(unsigned int kel=0; kel<AllElectrons.size(); kel++){
        for(unsigned int jk=0; jk < AK4PuppiJets->size(); jk++){
          bool isElectronFound = false;
          double this_dR_AKel = deltaR(AK4PuppiJets->at(jk).eta(), AK4PuppiJets->at(jk).phi(), AllElectrons[kel].eta(), AllElectrons[kel].phi());
            if(this_dR_AKel<0.4 && !isElectronFound){
              AK4lep_pt.push_back(AllElectrons[kel].pt());
              AK4lep_eta.push_back(AllElectrons[kel].eta());
              AK4lep_phi.push_back(AllElectrons[kel].phi());
              AK4lep_mass.push_back(AllElectrons[kel].mass());
              AK4lep_id.push_back(AllElectrons[kel].pdgId());
              isElectronFound = true;//stop jet cicle if the muon is asscociated to one AK4 jet
            }
          }
        }

      //Zbb event selection
      /*if(leadingAK8_pt_idx>-1&& subleadingAK8_pt_idx>-1){
	if(AK8PuppiJets->at(leadingAK8_pt_idx).pt()>450 && abs(AK8PuppiJets->at(leadingAK8_pt_idx).eta())<2.4 && AK8PuppiJets->at(leadingAK8_pt_idx).userFloat("ak8PFJetsPuppiSoftDropMass")>80 && AK8PuppiJets->at(leadingAK8_pt_idx).userFloat("ak8PFJetsPuppiSoftDropMass")<110 && AK8PuppiJets->at(subleadingAK8_pt_idx).pt()>200 && abs(AK8PuppiJets->at(subleadingAK8_pt_idx).eta())<2.4){
                if(Nmu==0 && Ne==0){passedZqqSelection=true;}
             }
      }*/
     //end Zbb event selection


     //hlt jets
     //std::cout<<"hltPFJetForBtag size: "<< hltjets->size()<<std::endl;
     //std::cout<<"pfJetTagCollection: "<<pfJetTagCollection->size()<<std::endl;
     /*for(unsigned int ijet=0; ijet<hltjetsForBTag->size(); ijet++){
     //std::cout<<"index jet: "<<ijet<<std::endl;
     //std::cout<<"jet pt: "<<hltjets->at(ijet).pt()<<std::endl;
     hltjetForBTag_pt.push_back(hltjetsForBTag->at(ijet).pt());
     hltjetForBTag_eta.push_back(hltjetsForBTag->at(ijet).eta());
     hltjetForBTag_phi.push_back(hltjetsForBTag->at(ijet).phi());
     hltjetForBTag_mass.push_back(hltjetsForBTag->at(ijet).mass());
     //hltParticleNetONNXJetTags_probb.push_back(hltjets->at(ijet).bDiscriminator("hltParticleNetONNXJetTags:probb"));
     //hltParticleNetONNXJetTags_probc.push_back(hltjets->at(ijet).bDiscriminator("hltParticleNetONNXJetTags:probc"));
     //hltParticleNetONNXJetTags_probuds.push_back(hltjets->at(ijet).bDiscriminator("hltParticleNetONNXJetTags:probuds"));
     //hltParticleNetONNXJetTags_probtauh.push_back(hltjets->at(ijet).bDiscriminator("hltParticleNetONNXJetTags:probtauh"));
		
     float tagValue_b = -20;
     float tagValue_c = -20;
     float tagValue_uds = -20;
     float tagValue_g = -20;
     float tagValue_tauh = -20;
     float minDR2_b = 0.01;
     float minDR2_c = 0.01;
     float minDR2_uds = 0.01;
     float minDR2_g = 0.01;
     float minDR2_tauh = 0.01;
      
     int index_tag=0;
     //std::cout<<"pfJetTagCollection: "<<pfJetTagCollection->size()<<std::endl;
	
     for (auto const &tag : *pfJetTagCollectionParticleNetprobc) {
        float dR2 = reco::deltaR2(hltjetsForBTag->at(ijet), *(tag.first));
        //std::cout<<"tag "<<index_tag<<"   deltaR= "<<dR2<<std::endl;
        if (dR2 < minDR2_c) {
          minDR2_c = dR2;
          tagValue_c = tag.second;
        }
        index_tag++;
     }
     for (auto const &tag : *pfJetTagCollectionParticleNetprobb) {
       float dR2 = reco::deltaR2(hltjetsForBTag->at(ijet), *(tag.first));
       //std::cout<<"tag "<<index_tag<<"   deltaR= "<<dR2<<std::endl;
       if (dR2 < minDR2_b) {
         minDR2_b = dR2;
         tagValue_b = tag.second;
       }
       //index_tag++;
      }
      for (auto const &tag : *pfJetTagCollectionParticleNetprobuds) {
        float dR2 = reco::deltaR2(hltjetsForBTag->at(ijet), *(tag.first));
        //std::cout<<"tag "<<index_tag<<"   deltaR= "<<dR2<<std::endl;
        if (dR2 < minDR2_uds) {
          minDR2_uds = dR2;
          tagValue_uds = tag.second;
        }
        //index_tag++;
       }
      for (auto const &tag : *pfJetTagCollectionParticleNetprobg) {
        float dR2 = reco::deltaR2(hltjetsForBTag->at(ijet), *(tag.first));
        //std::cout<<"tag "<<index_tag<<"   deltaR= "<<dR2<<std::endl;
        if (dR2 < minDR2_g) {
          minDR2_g = dR2;
          tagValue_g = tag.second;
        }
        //index_tag++;
      }
      for (auto const &tag : *pfJetTagCollectionParticleNetprobtauh) {
        float dR2 = reco::deltaR2(hltjetsForBTag->at(ijet), *(tag.first));
        //std::cout<<"tag "<<index_tag<<"   deltaR= "<<dR2<<std::endl;
          if (dR2 < minDR2_tauh) {
            minDR2_tauh = dR2;
            tagValue_tauh = tag.second;
          }
	      //index_tag++;
      }
      hltParticleNetONNXJetTags_probc.push_back(tagValue_c);	
      hltParticleNetONNXJetTags_probb.push_back(tagValue_b);	
      hltParticleNetONNXJetTags_probuds.push_back(tagValue_uds);	
      hltParticleNetONNXJetTags_probg.push_back(tagValue_g);	
      hltParticleNetONNXJetTags_probtauh.push_back(tagValue_tauh);	
    } */ 

    //std::cout<<"hltPFJetForBtag size: "<< hltjets->size()<<std::endl;

    /*for (auto const &hltjet : *hltjets) {
      std::cout<<"jet pt "<<hltjet.pt()<<std::endl;
      hltjet_pt.push_back(hltjet.pt());
      float tagValue = -20;
      float minDR2 = 0.01;
      
      for (auto const &tag : *pfJetTagCollection) {
        float dR2 = reco::deltaR2(hltjet, *(tag.first));
          if (dR2 < minDR2) {
            minDR2 = dR2;
            tagValue = tag.second;
          }
      }
      hltParticleNetONNXJetTags_probc.push_back(tagValue);	
      
				
    }*/
}


void HccAna::setGENVariables(edm::Handle<reco::GenParticleCollection> prunedgenParticles,
                                 edm::Handle<edm::View<pat::PackedGenParticle> > packedgenParticles,
                                 edm::Handle<edm::View<reco::GenJet> > genJets)
{
  reco::GenParticleCollection::const_iterator genPart;
  bool first_quarkgen=false;
  if(isHcc){
    for(genPart = prunedgenParticles->begin(); genPart != prunedgenParticles->end(); genPart++) {
      if(abs(genPart->pdgId())==1 || abs(genPart->pdgId())==2 || abs(genPart->pdgId())==3 || abs(genPart->pdgId())==4 || abs(genPart->pdgId())==5 || abs(genPart->pdgId())==6  || abs(genPart->pdgId())==7  || abs(genPart->pdgId())==23  || abs(genPart->pdgId())==24 || abs(genPart->pdgId())==25){
        //cout<<"pdg: "<< abs(genPart->pdgId()) <<"  pT: "<<genPart->pt() <<"   eta: "<<genPart->eta() <<"   phi: "<<genPart->phi() <<endl;
        const reco::Candidate * mom = genPart->mother(0);
        //cout<<"mother: "<< mom->pdgId()<<"   pT: "<< mom->pt() <<"   eta: "<< mom->eta() <<"   phi: "<< mom->phi() <<endl;
        bool Higgs_daughter=false;
        int n = genPart->numberOfDaughters();
        if(mom->pdgId()==2212 && first_quarkgen==false){
        for(int j_d = 0; j_d < n; ++ j_d) {
          const reco::Candidate * d = genPart->daughter( j_d );
          if((d->pdgId())==25){
            Higgs_daughter=true;
          }

          if(Higgs_daughter==true && (abs(d->pdgId())==1 || abs(d->pdgId())==2 || abs(d->pdgId())==3  || abs(d->pdgId())==4 || abs(d->pdgId())==5 || abs(d->pdgId())==6 || abs(d->pdgId())==7)){
            quark_pt.push_back(d->pt());
            quark_eta.push_back(d->eta());
            quark_phi.push_back(d->phi());
            quark_flavour.push_back(d->pdgId());
            quark_VBF.push_back(true);
            first_quarkgen=true;
          }

        }
      }

      if(( abs(genPart->pdgId())==4 || abs(genPart->pdgId())==5) && (mom->pdgId())==25){
        quark_pt.push_back(genPart->pt());
        quark_eta.push_back(genPart->eta());
        quark_phi.push_back(genPart->phi());
        quark_flavour.push_back(genPart->pdgId());
        quark_VBF.push_back(false);
       }

     }

   }

    /*if( abs(genPart->pdgId())==1 || abs(genPart->pdgId())==2 || abs(genPart->pdgId())==3 || abs(genPart->pdgId())==4 || abs(genPart->pdgId())==5 || abs(genPart->pdgId())==17 || abs(genPart->pdgId())==21   ){
    quark_pt.push_back(genPart->pt());
    quark_eta.push_back(genPart->eta());
    quark_phi.push_back(genPart->phi());
    quark_flavour.push_back(genPart->pdgId());
    }*/

  } // end isHcc condition

  if(isZcc){
    for(genPart = prunedgenParticles->begin(); genPart != prunedgenParticles->end(); genPart++) {
      if(genPart->pdgId()==23 && genPart->numberOfDaughters()==2){
        //std::cout<<"found Z candidate"<<std::endl;
        const reco::Candidate * da1 = genPart->daughter(0);
        const reco::Candidate * da2 = genPart->daughter(1);
        //std::cout<<"daughters: "<<da1->pdgId()<<"   "<<da2->pdgId()<<std::endl;
        if(fabs(da1->pdgId())==fabs(da2->pdgId()) && fabs(da1->pdgId())==4){
          //c1
          //std::cout<<"Zcc"<<endl;
          quark_pt.push_back(da1->pt());
          quark_eta.push_back(da1->eta());
          quark_phi.push_back(da1->phi());
          quark_flavour.push_back(da1->pdgId());
          //c2
          quark_pt.push_back(da2->pt());
          quark_eta.push_back(da2->eta());
          quark_phi.push_back(da2->phi());
          quark_flavour.push_back(da2->pdgId());
        }
      }
    }
  }
  
  if(isZqq){
    for(genPart = prunedgenParticles->begin(); genPart != prunedgenParticles->end(); genPart++) {
     if(abs(genPart->pdgId())==1 || abs(genPart->pdgId())==2 || abs(genPart->pdgId())==3 || abs(genPart->pdgId())==4 || abs(genPart->pdgId())==5 || abs(genPart->pdgId())==6 ){
       const reco::Candidate * Mom = genPart->mother(0);
       if( Mom->pdgId()==23){
            quark_pt.push_back(genPart->pt());
            quark_eta.push_back(genPart->eta());
            quark_phi.push_back(genPart->phi());
            quark_flavour.push_back(genPart->pdgId());
            // if the procedure is correct i must have only one Z boson,
            if(Z_pt != -1000 && Z_eta != -1000 && Z_phi != -1000 && Z_mass != -1000 && Z_pt != Mom->pt() && Z_eta != Mom->eta() && Z_phi != Mom->phi() && Z_mass != Mom->mass()){cerr << "Error in the event you are filling two different Z bosons!!";}
            Z_pt = Mom->pt();
            Z_eta = Mom->eta();
            Z_phi = Mom->phi();
            Z_mass = Mom->mass();
       } 
     }
    }
  }

	edm::View<reco::GenJet>::const_iterator genjet;

  for(genjet = genJets->begin(); genjet != genJets->end(); genjet++) {

	GENjet_pt.push_back(genjet->pt());
	GENjet_eta.push_back(genjet->eta());
	GENjet_phi.push_back(genjet->phi());
	GENjet_mass.push_back(genjet->mass());
	} //loop over gen jets


}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HccAna::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HccAna);

//  LocalWords:  ecalDriven
