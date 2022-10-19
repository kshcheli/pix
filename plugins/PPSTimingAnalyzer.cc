// -*- C++ -*-
//
// Package:    PPSTimingAnalyzer/PPSTimingAnalyzer
// Class:      PPSTimingAnalyzer
//
/**\class PPSTimingAnalyzer PPSTimingAnalyzer.cc PPSTimingAnalyzer/PPSTimingAnalyzer/plugins/PPSTimingAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Simple ntuple maker for PPS timing detector studies. Runs on AOD, and writes flat ntuples with arrays of RecHit, local track, and proton object (if available) variables. 
// Also saves information on central CMS vertices for low-pileup studies
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/Registry.h"

#include "FWCore/Framework/interface/EventSetup.h"

// HLT information                                   
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

// Central                                                                                                                                                                                   
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

// RP                                                                                                                                                                                        
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"

// Timing                                                                                                                                                                                    
#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSDiamondRecHit.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

// Protons                                                                                                                                                                                   
#include "DataFormats/ProtonReco/interface/ForwardProton.h"
#include "DataFormats/ProtonReco/interface/ForwardProtonFwd.h"

#include "CondFormats/RunInfo/interface/LHCInfo.h"
#include "CondFormats/DataRecord/interface/LHCInfoRcd.h"
#include "DataFormats/CTPPSReco/interface/CTPPSPixelRecHit.h"

#include "CondFormats/PPSObjects/interface/CTPPSPixelIndices.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"
#include "DataFormats/CTPPSDigi/interface/CTPPSPixelDigi.h"
#include "DataFormats/CTPPSDigi/interface/CTPPSPixelDataError.h"
#include "DataFormats/CTPPSReco/interface/CTPPSPixelCluster.h"
#include "DataFormats/CTPPSReco/interface/CTPPSPixelLocalTrack.h"
#include "DataFormats/CTPPSReco/interface/CTPPSPixelRecHit.h"




#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TH1.h>
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"


#include <TH2.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TSpline.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TF1.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <set>

#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TLorentzVector.h>


using namespace std;
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class PPSTimingAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit PPSTimingAnalyzer(const edm::ParameterSet&);
  ~PPSTimingAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  // Input tags                                                                                                                                                                               

  std::string outputFile_;
  unsigned int maxVertices_;

  edm::EDGetTokenT< edm::DetSetVector<CTPPSDiamondRecHit> > tokenDiamondHit_;
  edm::EDGetTokenT<std::vector<CTPPSLocalTrackLite> > pps_tracklite_token_;


  edm::EDGetTokenT< edm::View<reco::Vertex> > verticesToken_;
  edm::EDGetTokenT<reco::TrackCollection>  tracksToken_;
  edm::EDGetTokenT<edm::DetSetVector<CTPPSPixelRecHit>> tokenRecHit_;

  edm::EDGetTokenT<std::vector<reco::ForwardProton> > recoProtonsSingleRPToken_;
  edm::EDGetTokenT<std::vector<reco::ForwardProton> > recoProtonsMultiRPToken_;
  

  edm::ESGetToken<LHCInfo, LHCInfoRcd> lhcInfoToken_;

  std::string lhcInfoLabel_;

  // Tree contents                                                                                                                                                           

  // Run+event quantities                                                                                                                                                                   
  Int_t BX, Run, LumiSection, EventNum;
  Float_t CrossingAngle;

  Double_t TimingRecHitT[100], TimingRecHitX[100], TimingRecHitY[100], TimingRecHitToT[100];
  unsigned int TimingRecHitChannel[100], TimingRecHitArm[100], TimingRecHitPlane[100], TimingRecHitOOTIndex[100], TimingRecHitMultiHit[100];

  Int_t nVertices, nArmsTiming, nTracksTiming, nRecHitsTiming;

  Int_t nLiteTracks;
  Float_t TrackLiteX[100], TrackLiteY[100], TrackLiteZ[100], TrackLiteTime[100], TrackLiteTimeUnc[100];
  Int_t TrackLiteRPID[100], TrackLiteArm[100];

  Double_t PrimVertexZ[100], PrimVertexY[100], PrimVertexX[100];
  Int_t PrimVertexIsBS[100], PrimVertexNtracks[100];

  Int_t nProtons;
  Float_t ProtonXi[100];
  Float_t ProtonThY[100];
  Float_t ProtonThX[100];
  Float_t Protont[100];
  Int_t ProtonIsMultiRP[100];
  Int_t ProtonRPID[100];
  Int_t ProtonArm[100];
  Float_t ProtonTime[100];



  static constexpr int NArms = 2;
  static constexpr int NStationMAX = 3;  // in an arm
  static constexpr int NRPotsMAX = 6;    // per station
  static constexpr int NplaneMAX = 6;    // per RPot
  static constexpr int NROCsMAX = 6;     // per plane
  static constexpr int RPn_first = 3, RPn_last = 4;
  static constexpr int ADCMax = 256;
  static constexpr int StationIDMAX = 4;  // possible range of ID
  static constexpr int RPotsIDMAX = 8;    // possible range of ID
  static constexpr int NLocalTracksMAX = 20;
  static constexpr int hitMultMAX = 50;   // tuned
  static constexpr int ClusMultMAX = 10;  // tuned
  static constexpr int ClusterSizeMax = 9;
  static constexpr int mapXbins = 200;
  static constexpr int mapYbins = 240;
  static constexpr float mapYmin = -16.;
  static constexpr float mapYmax = 8.;
  const float mapXmin = 0. * TMath::Cos(18.4 / 180. * TMath::Pi());
  const float mapXmax = 30. * TMath::Cos(18.4 / 180. * TMath::Pi());
  int TrackFitDimension = 4;
  static constexpr int NRPotBinsInStation = RPn_last - RPn_first;
  static constexpr int NPlaneBins = NplaneMAX * NRPotBinsInStation;
  static constexpr int RPotsTotalNumber = NArms * NStationMAX * NRPotsMAX;

  TH1D *hBX[RPotsTotalNumber][NplaneMAX];
  
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
TFile* file;
TTree* tree;

//
// constructors and destructor
//
PPSTimingAnalyzer::PPSTimingAnalyzer(const edm::ParameterSet& iConfig)
//    : tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))) {
  : tokenDiamondHit_  ( consumes< edm::DetSetVector<CTPPSDiamondRecHit> >    ( iConfig.getParameter<edm::InputTag>( "tagDiamondRecHits" ) ) ),
    pps_tracklite_token_ ( consumes<std::vector<CTPPSLocalTrackLite>>(iConfig.getParameter<edm::InputTag>("tagTrackLites") ) ),
    verticesToken_    ( consumes< edm::View<reco::Vertex> >( iConfig.getParameter<edm::InputTag>( "verticesTag" ) ) ),
    tracksToken_    ( consumes< reco::TrackCollection >( iConfig.getParameter<edm::InputTag>( "tracksTag" ) ) ),
    tokenRecHit_    ( consumes< edm::DetSetVector<CTPPSPixelRecHit>>(iConfig.getParameter<edm::InputTag>( "tagRPixRecHit" ) ) ),
    recoProtonsSingleRPToken_   ( consumes<std::vector<reco::ForwardProton> >      ( iConfig.getParameter<edm::InputTag>( "ppsRecoProtonSingleRPTag" ) ) ),
    recoProtonsMultiRPToken_   ( consumes<std::vector<reco::ForwardProton> >      ( iConfig.getParameter<edm::InputTag>( "ppsRecoProtonMultiRPTag" ) ) ) {


  //now do what ever initialization is needed
  outputFile_ = iConfig.getUntrackedParameter<std::string>("outfilename", "output.root");
  maxVertices_ = iConfig.getParameter<uint>("maxVertices");
  file = new TFile(outputFile_.c_str(), "recreate");
  file->cd();
  // tree definition                                                                                                                                                         
  tree = new TTree("ntp1", "ntp1");



//////
//////
//////
        
//cout << "Inititalization"<<endl;
TString nom;
        for (int rrp = 0; rrp < RPotsTotalNumber; rrp++) {
        for (int ppp = 0; ppp < NplaneMAX; ppp++) {

  nom.Form("RP %d, plane %d, Nrechit vs BX", rrp,ppp);

//cout << "Inititalization2"<<endl;

//            hBX[rrp][ppp]->SetTitle(nom);

            hBX[rrp][ppp]= new TH1D(nom,nom ,4002, -1.5, 4000. + 0.5); //added        
//cout << "Inititalization3"<<endl;

          }
        }  // end of for(int p=0; p<NplaneMAX;..


///////
//////
//////

}

PPSTimingAnalyzer::~PPSTimingAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

for(int i=0;i<RPotsTotalNumber; i++){
 for(int j=0;j<NplaneMAX;j++){

  if(hBX[i][j]->GetEntries()!=0){

  hBX[i][j]->Write();
  }
 }
}


  file->Write();
  file->Close();

//cout << "SOMETHING3"<<endl;
}

//
// member functions
//

// ------------ method called for each event  ------------
void PPSTimingAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  // Run and BX information                                                                                                                          
  BX = iEvent.bunchCrossing();
  Run = iEvent.id().run();
  LumiSection = iEvent.luminosityBlock();
  EventNum = iEvent.id().event();



//cout << "SOMETHING1"<<endl;

  // PPS timing detector RecHits (1/channel)                                                                                                         
//  edm::Handle< edm::DetSetVector<CTPPSDiamondRecHit> > diamondRecHits;
 // iEvent.getByToken( tokenDiamondHit_, diamondRecHits );

  // Vertices in the central CMS tracker                
 // edm::Handle< edm::View<reco::Vertex> > vertices;
 // iEvent.getByToken( verticesToken_, vertices );

  // First filter on number of vertices for Low-PU data analysis
//  if(vertices->size() > maxVertices_)
//cout << "Run selection"<<endl;

///if(Run != 355933) return;
//cout << "After run sel"<<endl;
/*
//    return;

  // Tracks in the central CMS tracker                                                                                                               
  edm::Handle<reco::TrackCollection> tks;
  iEvent.getByToken(tracksToken_, tks);

  // Initializing counters                                                                                                                           
  nProtons = 0;
  nTracksTiming = 0;
  nRecHitsTiming = 0;
  nLiteTracks = 0;
  nVertices = 0;

  // Initializing arrays                                                                                                                             
  for(int i = 0; i < 100; i++)
    {
      PrimVertexZ[i] = -999.;
      PrimVertexX[i] = -999.;
      PrimVertexY[i] = -999.;
      PrimVertexIsBS[i] = -999;
      PrimVertexNtracks[i] = -999;

      TimingRecHitT[i] = -999;
      TimingRecHitX[i] = -999;
      TimingRecHitY[i] = -999;
      TimingRecHitOOTIndex[i] = -999;
      TimingRecHitMultiHit[i] = -999;
      TimingRecHitToT[i] = 0;
      TimingRecHitChannel[i] = -999;
      TimingRecHitArm[i] = -999;
      TimingRecHitPlane[i] = -999;

      ProtonXi[i] = -999;
      ProtonThY[i] = -999;
      ProtonThX[i] = -999;
      Protont[i] = -999;
      ProtonIsMultiRP[i] = -999;
      ProtonRPID[i] = -999;
      ProtonArm[i] = -999;
      ProtonTime[i] = -999.;
      TrackLiteX[i] = 0;
      TrackLiteY[i] = 0;
      TrackLiteZ[i] = 0;
      TrackLiteArm[i] = 0;
      TrackLiteTime[i] = 0;
      TrackLiteTimeUnc[i] = 0;
      TrackLiteRPID[i] = 0;
    }

  // RecHits - timing 
  for ( const auto& rechits_ds : *diamondRecHits ) {
    const CTPPSDiamondDetId detidforrh( rechits_ds.detId() );
    for ( const auto& rechit : rechits_ds ) {

      TimingRecHitArm[nRecHitsTiming] = detidforrh.arm();
      TimingRecHitChannel[nRecHitsTiming] = detidforrh.channel();
      TimingRecHitPlane[nRecHitsTiming] = detidforrh.plane();
      TimingRecHitT[nRecHitsTiming] = rechit.time();
      TimingRecHitX[nRecHitsTiming] = rechit.x();
      TimingRecHitY[nRecHitsTiming] = rechit.y();
      TimingRecHitToT[nRecHitsTiming] = rechit.toT();
      TimingRecHitOOTIndex[nRecHitsTiming] = rechit.ootIndex();                                                                       
      TimingRecHitMultiHit[nRecHitsTiming] = rechit.multipleHits();

      nRecHitsTiming++;
cout << "timingrechit"<<endl;
    }
  }
//cout << "SOMETHING"<<endl;
  // Proton lite tracks 
  edm::Handle<std::vector<CTPPSLocalTrackLite> > ppsTracksLite;
  iEvent.getByToken( pps_tracklite_token_, ppsTracksLite );

  for ( const auto& trklite : *ppsTracksLite )
    {
      const CTPPSDetId detid( trklite.rpId() );

      // transform the raw, 32-bit unsigned integer detId into the TOTEM "decimal" notation                                                            
      const unsigned short raw_id = 100*detid.arm()+10*detid.station()+detid.rp();

      TrackLiteX[nLiteTracks] = trklite.x();
      TrackLiteY[nLiteTracks] = trklite.y();
      TrackLiteTime[nLiteTracks] = trklite.time();
      TrackLiteTimeUnc[nLiteTracks] = trklite.timeUnc();
      TrackLiteRPID[nLiteTracks] = raw_id;

      if(raw_id == 3)
	{
	  TrackLiteZ[nLiteTracks] = 210;
	  TrackLiteArm[nLiteTracks] = 0;
	}
      if(raw_id == 23)
        {
          TrackLiteZ[nLiteTracks] = 220;
          TrackLiteArm[nLiteTracks] = 0;
	}
      if(raw_id == 16)
        {
          TrackLiteZ[nLiteTracks] = 216;
          TrackLiteArm[nLiteTracks] = 0;
        }
      if(raw_id == 103)
        {
          TrackLiteZ[nLiteTracks] = 210;
          TrackLiteArm[nLiteTracks] = 1;
        }
      if(raw_id == 123)
        {
          TrackLiteZ[nLiteTracks] = 220;
          TrackLiteArm[nLiteTracks] = 1;
        }
      if(raw_id == 116)
        {
          TrackLiteZ[nLiteTracks] = 216;
          TrackLiteArm[nLiteTracks] = 1;
        }


      nLiteTracks++;
    }

  // Full Reco protons 
  edm::Handle<std::vector<reco::ForwardProton>> recoMultiRPProtons;
  iEvent.getByToken(recoProtonsMultiRPToken_, recoMultiRPProtons);
  edm::Handle<std::vector<reco::ForwardProton>> recoSingleRPProtons;
  iEvent.getByToken(recoProtonsSingleRPToken_, recoSingleRPProtons);

  int ismultirp = -999;
  unsigned int decRPId = -999;
  unsigned int armId = -999;
  float th_y = -999;
  float th_x = -999;
  float t = -999;
  float xi = -999.;
  float protontime = -999.;

  // Proton object with Single-RP algorithm                                                                                                          
  for (const auto & proton : *recoSingleRPProtons)
    {
      if (proton.validFit())
        {
          th_y = proton.thetaY();
          th_x = proton.thetaX();
          xi = proton.xi();
          //      t = proton.t();                                                                                                                      
          protontime = proton.time();

          CTPPSDetId rpId((*proton.contributingLocalTracks().begin())->rpId());
          decRPId = rpId.arm()*100 + rpId.station()*10 + rpId.rp();

          ismultirp = 0;

          ProtonXi[nProtons] = xi;
          ProtonThX[nProtons] = th_x;
          ProtonThY[nProtons] = th_y;
          Protont[nProtons] = t;
          ProtonIsMultiRP[nProtons] = ismultirp;
          ProtonRPID[nProtons] = decRPId;
          ProtonArm[nProtons] = armId;
          ProtonTime[nProtons] = protontime;
          nProtons++;
        }
    }

  // Proton object with Multi-RP algorithm                                                                                                          
  for (const auto & proton : *recoMultiRPProtons)
    {
      if (proton.validFit())
        {
          th_y = proton.thetaY();
          th_x = proton.thetaX();
          xi = proton.xi();
          //      t = proton.t();                                                                                                                    
          protontime = proton.time();

          CTPPSDetId rpId((*proton.contributingLocalTracks().begin())->rpId());
          armId = rpId.arm();

          ismultirp = 1;

          ProtonXi[nProtons] = xi;
          ProtonThX[nProtons] = th_x;
          ProtonThY[nProtons] = th_y;
          Protont[nProtons] = t;
          ProtonIsMultiRP[nProtons] = ismultirp;
          ProtonRPID[nProtons] = decRPId;
          ProtonArm[nProtons] = armId;
          ProtonTime[nProtons] = protontime;
          nProtons++;
        }
    }

  // Primary vertices in central CMS tracker
  for ( const auto& vtx : *vertices ) {
    PrimVertexZ[nVertices] = vtx.z();
    PrimVertexX[nVertices] = vtx.x();
    PrimVertexY[nVertices] = vtx.y();
    PrimVertexNtracks[nVertices] = vtx.nTracks(0);

    if(vtx.isFake() == 1)
      PrimVertexIsBS[nVertices] = 1;
    else
      PrimVertexIsBS[nVertices] = 0;
    nVertices++;
  }
*/
///////////////////////////////////////////////
//cout << "bef rechit "<<endl;
 Handle<edm::DetSetVector<CTPPSPixelRecHit>> pixRecHits;
  iEvent.getByToken(tokenRecHit_, pixRecHits);

  for (auto &recHitDs : *pixRecHits) {    
    CTPPSPixelDetId recHitId = CTPPSPixelDetId(recHitDs.id);

  cout << " ~~~> Found recHit! "<< recHitId.rp()<<" " << recHitId.plane()<<endl;
    hBX[recHitId.rp()][recHitId.plane()]->Fill(BX);

 
  }


	edm::DetSetVector<CTPPSPixelRecHit>::const_iterator rhDSViter = pixRecHits->begin();
for (; rhDSViter != pixRecHits->end(); rhDSViter++) {
edm::DetSet<CTPPSPixelRecHit>::const_iterator begin = (*rhDSViter).begin();
				edm::DetSet<CTPPSPixelRecHit>::const_iterator end = (*rhDSViter).end();
				for (edm::DetSet<CTPPSPixelRecHit>::const_iterator rh = begin; rh != end; rh++) {

    CTPPSPixelDetId recHitIdt = CTPPSPixelDetId(rhDSViter->detId() );
 cout  <<"SECOND ~~~> Found recHit! "<< recHitIdt.rp()<<" " << recHitIdt.plane()<<endl;
 
 
  }}


//cout << "after rechit "<<endl;

//////////////////////////////////////////////

  //tree->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void PPSTimingAnalyzer::beginJob() {
  // Booking the ntuple                                                                                                                                

  tree->Branch("Run", &Run, "Run/I");
  tree->Branch("LumiSection", &LumiSection, "LumiSection/I");
  tree->Branch("BX", &BX, "BX/I");
  tree->Branch("EventNum", &EventNum, "EventNum/I");
  tree->Branch("CrossingAngle", &CrossingAngle, "CrossingAngle/F");

  tree->Branch("nRecHitsTiming", &nRecHitsTiming, "nRecHitsTiming/I");
  tree->Branch("TimingRecHitArm",&TimingRecHitArm,"TimingRecHitArm[nRecHitsTiming]/I");
  tree->Branch("TimingRecHitChannel",&TimingRecHitChannel,"TimingRecHitChannel[nRecHitsTiming]/I");
  tree->Branch("TimingRecHitPlane",&TimingRecHitPlane,"TimingRecHitPlane[nRecHitsTiming]/I");
  tree->Branch("TimingRecHitT",&TimingRecHitT,"TimingRecHitT[nRecHitsTiming]/D");
  tree->Branch("TimingRecHitX",&TimingRecHitX,"TimingRecHitX[nRecHitsTiming]/D");
  tree->Branch("TimingRecHitY",&TimingRecHitY,"TimingRecHitY[nRecHitsTiming]/D");
  tree->Branch("TimingRecHitToT",&TimingRecHitToT,"TimingRecHitToT[nRecHitsTiming]/D");
  tree->Branch("TimingRecHitOOTIndex",&TimingRecHitOOTIndex,"TimingRecHitOOTIndex[nRecHitsTiming]/I");                                      
  tree->Branch("TimingRecHitMultiHit",&TimingRecHitMultiHit,"TimingRecHitMultiHit[nRecHitsTiming]/I");                                      

  tree->Branch("nLiteTracks", &nLiteTracks, "nLiteTracks/I");
  tree->Branch("TrackLiteX", &TrackLiteX, "TrackLiteX[nLiteTracks]/F");
  tree->Branch("TrackLiteY", &TrackLiteY, "TrackLiteY[nLiteTracks]/F");
  tree->Branch("TrackLiteZ", &TrackLiteZ, "TrackLiteZ[nLiteTracks]/F");
  tree->Branch("TrackLiteTime", &TrackLiteTime, "TrackLiteTime[nLiteTracks]/F");
  tree->Branch("TrackLiteTimeUnc", &TrackLiteTimeUnc, "TrackLiteTimeUnc[nLiteTracks]/F");
  tree->Branch("TrackLiteRPID", &TrackLiteRPID, "TrackLiteRPID[nLiteTracks]/I");
  tree->Branch("TrackLiteArm", &TrackLiteArm, "TrackLiteArm[nLiteTracks]/I");

  tree->Branch("nVertices", &nVertices, "nVertices/I");
  tree->Branch("PrimVertexZ", &PrimVertexZ, "PrimVertexZ[nVertices]/D");
  tree->Branch("PrimVertexX", &PrimVertexX, "PrimVertexX[nVertices]/D");
  tree->Branch("PrimVertexY", &PrimVertexY, "PrimVertexY[nVertices]/D");
  tree->Branch("PrimVertexNtracks", &PrimVertexNtracks, "PrimVertexNtracks[nVertices]/I");
  tree->Branch("PrimVertexIsBS", &PrimVertexIsBS, "PrimVertexIsBS[nVertices]/I");

  tree->Branch("nProtons", &nProtons, "nProtons/I");
  tree->Branch("ProtonXi", &ProtonXi, "ProtonXi[nProtons]/F");
  tree->Branch("ProtonThX", &ProtonThX, "ProtonThX[nProtons]/F");
  tree->Branch("ProtonThY", &ProtonThY, "ProtonThY[nProtons]/F");
  tree->Branch("Protont", &Protont, "Protont[nProtons]/F");
  tree->Branch("ProtonIsMultiRP", &ProtonIsMultiRP, "ProtonIsMultiRP[nProtons]/I");
  tree->Branch("ProtonRPID", &ProtonRPID, "ProtonRPID[nProtons]/I");
  tree->Branch("ProtonArm", &ProtonArm, "ProtonArm[nProtons]/I");
  tree->Branch("ProtonTime", &ProtonTime, "ProtonTime[nProtons]/F");

}

// ------------ method called once each job just after ending the event loop  ------------
void PPSTimingAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void PPSTimingAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PPSTimingAnalyzer);
