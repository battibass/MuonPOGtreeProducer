//////////////////////////////////////
// Ntuplizer that fills muon_pog trees
//////////////////////////////////////

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Conversion.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/Common/interface/View.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"

#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/L1Trigger/interface/Muon.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenStatusFlags.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

#include "DataFormats/Scalers/interface/LumiScalers.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"

#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "MuonPOGtreeProducer/Tools/src/MuonPogTree.h"
#include "TTree.h"

#include <algorithm>
#include <iostream>

class MuonPogTreeProducer : public edm::EDAnalyzer 
{
public:

  MuonPogTreeProducer(const edm::ParameterSet &);
  
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void beginJob();
  virtual void endJob();
  
private:
  
  void fillGenInfo(const edm::Handle<std::vector<PileupSummaryInfo> > &,
		   const  edm::Handle<GenEventInfoProduct> &);

  void fillGenParticles(const edm::Handle<reco::GenParticleCollection> &);

  void fillHlt(const edm::Handle<edm::TriggerResults> &, 
	       const edm::Handle<trigger::TriggerEvent> &,
	       const edm::TriggerNames &);
  
  void fillPV(const edm::Handle<std::vector<reco::Vertex> > &);
  

  void fillDtSegments(const edm::Handle<DTRecSegment4DCollection> &);

  void fillCscSegments(const edm::Handle<CSCSegmentCollection> &);

  Int_t fillMuons(const edm::Handle<edm::View<reco::Muon> > &,
		  const edm::Handle<std::vector<reco::Vertex> > &,
		  const edm::Handle<reco::BeamSpot> &);

  void fillL1(const edm::Handle<l1t::MuonBxCollection> &);

    
  // returns false in case the match is for a RPC chamber
  bool getMuonChamberId(const DetId & id, muon_pog::MuonDetType & det, Int_t & r, Int_t & phi, Int_t & eta) const ;
  
  edm::EDGetTokenT<edm::TriggerResults> trigResultsToken_;
  edm::EDGetTokenT<trigger::TriggerEvent> trigSummaryToken_;

  std::string trigFilterCut_;
  std::string trigPathCut_;

  edm::EDGetTokenT<edm::View<reco::Muon> > muonToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex> > primaryVertexToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;

  edm::EDGetTokenT<DTRecSegment4DCollection> dtSegmentToken_;
  edm::EDGetTokenT<CSCSegmentCollection>     cscSegmentToken_;

  edm::EDGetTokenT<reco::PFMETCollection> pfMetToken_;
  edm::EDGetTokenT<reco::PFMETCollection> pfChMetToken_;
  edm::EDGetTokenT<reco::CaloMETCollection> caloMetToken_;

  edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileUpInfoToken_;
  edm::EDGetTokenT<GenEventInfoProduct> genInfoToken_;

  edm::EDGetTokenT<LumiScalersCollection> scalersToken_;
    
  edm::EDGetTokenT<l1t::MuonBxCollection> l1Token_;

  Float_t m_minMuPtCut;
  Int_t m_minNMuCut;

  muon_pog::Event event_;
  muon_pog::EventId eventId_;
  std::map<std::string,TTree*> tree_;
  
};


MuonPogTreeProducer::MuonPogTreeProducer( const edm::ParameterSet & cfg )
{

  // Input collections
  edm::InputTag tag = cfg.getUntrackedParameter<edm::InputTag>("TrigResultsTag", edm::InputTag("TriggerResults::HLT"));
  if (tag.label() != "none") trigResultsToken_ = consumes<edm::TriggerResults>(tag);

  tag = cfg.getUntrackedParameter<edm::InputTag>("TrigSummaryTag", edm::InputTag("hltTriggerSummaryAOD::HLT")); 
  if (tag.label() != "none") trigSummaryToken_ =consumes<trigger::TriggerEvent>(tag);

  trigFilterCut_ = cfg.getUntrackedParameter<std::string>("TrigFilterCut", std::string("all"));
  trigPathCut_ = cfg.getUntrackedParameter<std::string>("TrigPathCut", std::string("all"));

  tag = cfg.getUntrackedParameter<edm::InputTag>("DTSegmentTag", edm::InputTag("dt4DSegments"));
  if (tag.label() != "none") dtSegmentToken_ = consumes<DTRecSegment4DCollection>(tag);

  tag = cfg.getUntrackedParameter<edm::InputTag>("CSCSegmentTag", edm::InputTag("cscSegments"));
  if (tag.label() != "none") cscSegmentToken_ = consumes<CSCSegmentCollection>(tag);

  tag = cfg.getUntrackedParameter<edm::InputTag>("MuonTag", edm::InputTag("muons"));
  if (tag.label() != "none") muonToken_ = consumes<edm::View<reco::Muon> >(tag);

  tag = cfg.getUntrackedParameter<edm::InputTag>("MuonTag", edm::InputTag("muons"));
  if (tag.label() != "none") muonToken_ = consumes<edm::View<reco::Muon> >(tag);

  tag = cfg.getUntrackedParameter<edm::InputTag>("PrimaryVertexTag", edm::InputTag("offlinePrimaryVertices"));
  if (tag.label() != "none") primaryVertexToken_ = consumes<std::vector<reco::Vertex> >(tag);

  tag = cfg.getUntrackedParameter<edm::InputTag>("BeamSpotTag", edm::InputTag("offlineBeamSpot"));
  if (tag.label() != "none") beamSpotToken_ = consumes<reco::BeamSpot>(tag);

  tag = cfg.getUntrackedParameter<edm::InputTag>("PFMetTag", edm::InputTag("pfMet"));
  if (tag.label() != "none") pfMetToken_ = consumes<reco::PFMETCollection>(tag);

  tag = cfg.getUntrackedParameter<edm::InputTag>("PFChMetTag", edm::InputTag("pfChMet"));
  if (tag.label() != "none") pfChMetToken_ = consumes<reco::PFMETCollection>(tag);
 
  tag = cfg.getUntrackedParameter<edm::InputTag>("CaloMetTag", edm::InputTag("caloMet"));
  if (tag.label() != "none") caloMetToken_ = consumes<reco::CaloMETCollection>(tag); 

  tag = cfg.getUntrackedParameter<edm::InputTag>("GenTag", edm::InputTag("prunedGenParticles"));
  if (tag.label() != "none") genToken_ = consumes<reco::GenParticleCollection>(tag);

  tag = cfg.getUntrackedParameter<edm::InputTag>("PileUpInfoTag", edm::InputTag("pileupInfo"));
  if (tag.label() != "none") pileUpInfoToken_ = consumes<std::vector<PileupSummaryInfo> >(tag);

  tag = cfg.getUntrackedParameter<edm::InputTag>("GenInfoTag", edm::InputTag("generator"));
  if (tag.label() != "none") genInfoToken_ = consumes<GenEventInfoProduct>(tag);  

  tag = cfg.getUntrackedParameter<edm::InputTag>("ScalersTag", edm::InputTag("scalersRawToDigi"));
  if (tag.label() != "none") scalersToken_ = consumes<LumiScalersCollection>(tag);
    
  tag = cfg.getUntrackedParameter<edm::InputTag>("l1MuonsTag", edm::InputTag("gmtStage2Digis:Muon:"));
  if (tag.label() != "none") l1Token_ = consumes<l1t::MuonBxCollection>(tag);

  m_minMuPtCut = cfg.getUntrackedParameter<double>("MinMuPtCut", 0.);
  m_minNMuCut  = cfg.getUntrackedParameter<int>("MinNMuCut",  0.);

}


void MuonPogTreeProducer::beginJob() 
{
  
  edm::Service<TFileService> fs;
  tree_["muPogTree"] = fs->make<TTree>("MUONPOGTREE","Muon POG Tree");

  int splitBranches = 2;
  tree_["muPogTree"]->Branch("event",&event_,64000,splitBranches);
  tree_["muPogTree"]->Branch("eventId",&eventId_,64000,splitBranches);

}


void MuonPogTreeProducer::beginRun(const edm::Run & run, const edm::EventSetup & config )
{
  
}


void MuonPogTreeProducer::endJob() 
{

}


void MuonPogTreeProducer::analyze (const edm::Event & ev, const edm::EventSetup &)
{

  // Clearing branch variables
  // and setting default values
  event_.hlt.triggers.clear();
  event_.hlt.objects.clear();
  event_.l1muons.clear();

  event_.genParticles.clear();
  event_.genInfos.clear();
  event_.muons.clear();

  event_.dtSegments.clear();
  event_.cscSegments.clear();
  
  event_.mets.pfMet   = -999; 
  event_.mets.pfChMet = -999; 
  event_.mets.caloMet = -999; 

  for (unsigned int ix=0; ix<3; ++ix) {
    event_.primaryVertex[ix] = 0.;
    for (unsigned int iy=0; iy<3; ++iy) {
      event_.cov_primaryVertex[ix][iy] = 0.;
    }
  }
  event_.nVtx = -1;


  // Fill general information
  // run, luminosity block, event
  event_.runNumber = ev.id().run();
  event_.luminosityBlockNumber = ev.id().luminosityBlock();
  event_.eventNumber = ev.id().event();

  eventId_.runNumber = ev.id().run();
  eventId_.luminosityBlockNumber = ev.id().luminosityBlock();
  eventId_.eventNumber = ev.id().event();
    
  // Fill GEN pile up information
  if (!ev.isRealData()) 
    {
      if (!pileUpInfoToken_.isUninitialized() &&
	  !genInfoToken_.isUninitialized()) 
	{
	  edm::Handle<std::vector<PileupSummaryInfo> > puInfo;
	  edm::Handle<GenEventInfoProduct> genInfo;

	  if (ev.getByToken(pileUpInfoToken_, puInfo) &&
	      ev.getByToken(genInfoToken_, genInfo) ) 
	    fillGenInfo(puInfo,genInfo);
	  else 
	    edm::LogError("") << "[MuonPogTreeProducer]: Pile-Up Info collection does not exist !!!";
	}      
    }
  

  // Fill GEN particles information
  if (!ev.isRealData()) 
    {
      if (!genToken_.isUninitialized() ) 
	{ 
	  edm::Handle<reco::GenParticleCollection> genParticles;
	  if (ev.getByToken(genToken_, genParticles)) 
	    fillGenParticles(genParticles);
	  else 
	    edm::LogError("") << ">>> GEN collection does not exist !!!";
	}
    }

  if (ev.isRealData()) 
    {

      event_.bxId  = ev.bunchCrossing();
      event_.orbit = ev.orbitNumber();

      if (!scalersToken_.isUninitialized()) 
        { 
          edm::Handle<LumiScalersCollection> lumiScalers;
          if (ev.getByToken(scalersToken_, lumiScalers) && 
              lumiScalers->size() > 0 ) 
            event_.instLumi  = lumiScalers->begin()->instantLumi();
          else 
            edm::LogError("") << ">>> Scaler collection does not exist !!!";
        }
    }

  // Fill trigger information
  if (!trigResultsToken_.isUninitialized() &&
      !trigSummaryToken_.isUninitialized()) 
    {
      
      edm::Handle<edm::TriggerResults> triggerResults;
      edm::Handle<trigger::TriggerEvent> triggerEvent;
      
      if (ev.getByToken(trigResultsToken_, triggerResults) &&
	  ev.getByToken(trigSummaryToken_, triggerEvent)) 
	fillHlt(triggerResults, triggerEvent,ev.triggerNames(*triggerResults));
      else 
	edm::LogError("") << "[MuonPogTreeProducer]: Trigger collections do not exist !!!";
    }
  
  // Fill DT segment information
   if (!dtSegmentToken_.isUninitialized()) 
    {     

      edm::Handle<DTRecSegment4DCollection> dtSegments;

      if (ev.getByToken(dtSegmentToken_,dtSegments))
	fillDtSegments(dtSegments);
      else
	edm::LogError("") << "[MuonPogTreeProducer]: DT segments collection does not exist !!!";

    }

  // Fill CSC segment information
  if (!cscSegmentToken_.isUninitialized()) 
    {
      
      edm::Handle<CSCSegmentCollection> cscSegments;

      if (ev.getByToken(cscSegmentToken_,cscSegments))
	fillCscSegments(cscSegments);
      else
	edm::LogError("") << "[MuonPogTreeProducer]: CSC segments collection does not exist !!!";

    }

  // Fill vertex information
  edm::Handle<std::vector<reco::Vertex> > vertexes;

  if(!primaryVertexToken_.isUninitialized()) 
    {
      if (ev.getByToken(primaryVertexToken_, vertexes))
	fillPV(vertexes);
      else 
	edm::LogError("") << "[MuonPogTreeProducer]: Vertex collection does not exist !!!";
    }

  // Get beam spot for muons
  edm::Handle<reco::BeamSpot> beamSpot;
  if (!beamSpotToken_.isUninitialized() ) 
    { 
      if (!ev.getByToken(beamSpotToken_, beamSpot)) 
	edm::LogError("") << "[MuonPogTreeProducer]: Beam spot collection not found !!!";
    }

  // Fill (raw) MET information: PF, PF charged, Calo    
  edm::Handle<reco::PFMETCollection> pfMet; 
  if(!pfMetToken_.isUninitialized()) 
    { 
      if (!ev.getByToken(pfMetToken_, pfMet)) 
	edm::LogError("") << "[MuonPogTreeProducer] PFMet collection does not exist !!!"; 
      else { 
	const reco::PFMET &iPfMet = (*pfMet)[0]; 
	event_.mets.pfMet = iPfMet.et(); 
      } 
    } 

  edm::Handle<reco::PFMETCollection> pfChMet; 
  if(!pfChMetToken_.isUninitialized()) 
    { 
      if (!ev.getByToken(pfChMetToken_, pfChMet)) 
	edm::LogError("") << "[MuonPogTreeProducer] PFChMet collection does not exist !!!"; 
      else { 
	const reco::PFMET &iPfChMet = (*pfChMet)[0]; 
	event_.mets.pfChMet = iPfChMet.et(); 
      } 
    } 

  edm::Handle<reco::CaloMETCollection> caloMet; 
  if(!caloMetToken_.isUninitialized()) 
    { 
      if (!ev.getByToken(caloMetToken_, caloMet)) 
	edm::LogError("") << "[MuonPogTreeProducer] CaloMet collection does not exist !!!"; 
      else { 
	const reco::CaloMET &iCaloMet = (*caloMet)[0]; 
	event_.mets.caloMet = iCaloMet.et(); 
      } 
    } 

  // Get muons  
  edm::Handle<edm::View<reco::Muon> > muons;
  if (!muonToken_.isUninitialized() ) 
    { 
      if (!ev.getByToken(muonToken_, muons)) 
	edm::LogError("") << "[MuonPogTreeProducer] Muon collection does not exist !!!";
    }
  

  Int_t nGoodMuons = 0;
  eventId_.maxPTs.clear();
  // Fill muon information
  if (muons.isValid() && vertexes.isValid() && beamSpot.isValid()) 
    {
      nGoodMuons = fillMuons(muons,vertexes,beamSpot);
    }
  eventId_.nMuons = nGoodMuons;

    
  //Fill L1 informations
  edm::Handle<l1t::MuonBxCollection> l1s;
  if (!l1Token_.isUninitialized() )
    {
        if (!ev.getByToken(l1Token_, l1s))
	  edm::LogError("") << "[MuonPogTreeProducer] L1 muon bx collection does not exist !!!";
        else {
            fillL1(l1s);
        }
    }
    
  if (nGoodMuons >= m_minNMuCut)
  tree_["muPogTree"]->Fill();

}

void MuonPogTreeProducer::fillGenInfo(const edm::Handle<std::vector<PileupSummaryInfo> > & puInfo,
				      const edm::Handle<GenEventInfoProduct> & gen)
{

  muon_pog::GenInfo genInfo;
  
  genInfo.trueNumberOfInteractions     = -1.;
  genInfo.actualNumberOfInteractions   = -1.;
  genInfo.genWeight = gen->weight() ;
	
  std::vector<PileupSummaryInfo>::const_iterator puInfoIt  = puInfo->begin();
  std::vector<PileupSummaryInfo>::const_iterator puInfoEnd = puInfo->end();

  for(; puInfoIt != puInfoEnd; ++puInfoIt) 
    {
      int bx = puInfoIt->getBunchCrossing();
	  
      if(bx == 0) 
	{ 
	  genInfo.trueNumberOfInteractions   = puInfoIt->getTrueNumInteractions();
	  genInfo.actualNumberOfInteractions = puInfoIt->getPU_NumInteractions();
	  continue;
	}
    }
  
  event_.genInfos.push_back(genInfo);
  
}


void MuonPogTreeProducer::fillGenParticles(const edm::Handle<reco::GenParticleCollection> & genParticles)
{
  
  unsigned int gensize = genParticles->size();
  
  // Do not record the initial protons
  for (unsigned int i=0; i<gensize; ++i) 
    {

      const reco::GenParticle& part = genParticles->at(i);
    
      muon_pog::GenParticle gensel;
      gensel.pdgId = part.pdgId();
      gensel.status = part.status();
      gensel.energy = part.energy();
      gensel.pt = part.pt();
      gensel.eta = part.eta();
      gensel.phi = part.phi();
      gensel.vx = part.vx();
      gensel.vy = part.vy();
      gensel.vz = part.vz();

      // Full set of GenFlags
      gensel.flags.clear();
      reco::GenStatusFlags statusflags = part.statusFlags();
      if (statusflags.flags_.size() == 15)
	for (unsigned int flag = 0; flag < statusflags.flags_.size(); ++flag)
	  gensel.flags.push_back(statusflags.flags_[flag]);      
      
      gensel.mothers.clear();
      unsigned int nMothers = part.numberOfMothers();

      for (unsigned int iMother=0; iMother<nMothers; ++iMother) 
	{
	  gensel.mothers.push_back(part.motherRef(iMother)->pdgId());
	}

      // Protect agains bug in genParticles (missing mother => first proton)
      if (i>=2 && nMothers==0) gensel.mothers.push_back(0);
      
      event_.genParticles.push_back(gensel);
    }
  
}


void MuonPogTreeProducer::fillHlt(const edm::Handle<edm::TriggerResults> & triggerResults, 
				  const edm::Handle<trigger::TriggerEvent> & triggerEvent,
				  const edm::TriggerNames & triggerNames)
{    

  for (unsigned int iTrig=0; iTrig<triggerNames.size(); ++iTrig) 
    {
      
      if (triggerResults->accept(iTrig)) 
	{
	  std::string pathName = triggerNames.triggerName(iTrig);
	  if (trigPathCut_ == "all" || pathName.find(trigPathCut_) != std::string::npos)
	    event_.hlt.triggers.push_back(pathName);
	}
    }
      
  const trigger::size_type nFilters(triggerEvent->sizeFilters());

  for (trigger::size_type iFilter=0; iFilter!=nFilters; ++iFilter) 
    {
	
      std::string filterTag = triggerEvent->filterTag(iFilter).encode();
      
      if (trigFilterCut_ == "all" || filterTag.find(trigFilterCut_) != std::string::npos)
	{

	  trigger::Keys objectKeys = triggerEvent->filterKeys(iFilter);
	  const trigger::TriggerObjectCollection& triggerObjects(triggerEvent->getObjects());
	
	  for (trigger::size_type iKey=0; iKey<objectKeys.size(); ++iKey) 
	    {  
	      trigger::size_type objKey = objectKeys.at(iKey);
	      const trigger::TriggerObject& triggerObj(triggerObjects[objKey]);
	      
	      muon_pog::HLTObject hltObj;
	      
	      float trigObjPt = triggerObj.pt();
	      float trigObjEta = triggerObj.eta();
	      float trigObjPhi = triggerObj.phi();
	      
	      hltObj.filterTag = filterTag;
	      
	      hltObj.pt  = trigObjPt;
	      hltObj.eta = trigObjEta;
	      hltObj.phi = trigObjPhi;
	      
	      event_.hlt.objects.push_back(hltObj);
	      
	    }
	}
    }

}

void MuonPogTreeProducer::fillL1(const edm::Handle<l1t::MuonBxCollection> & l1MuonBxColl)
{

  for (int ibx = l1MuonBxColl->getFirstBX(); ibx <= l1MuonBxColl->getLastBX(); ++ibx) 
    {
      for (auto l1MuIt = l1MuonBxColl->begin(ibx); l1MuIt != l1MuonBxColl->end(ibx); ++l1MuIt)
	{

	  muon_pog::L1Muon l1part;
	  l1part.pt = l1MuIt->pt();
	  l1part.eta = l1MuIt->eta();
	  l1part.phi = l1MuIt->phi();
	  l1part.charge = l1MuIt->hwChargeValid() ? l1MuIt->charge() : 0;
	  
	  l1part.quality = l1MuIt->hwQual();
	  l1part.bx = ibx;
	  
	  l1part.tfIndex = l1MuIt->tfMuonIndex();
	  
	  event_.l1muons.push_back(l1part);
	  
	}
    }
}



void MuonPogTreeProducer::fillPV(const edm::Handle<std::vector<reco::Vertex> > & vertexes)
{
      
  int nVtx = 0;

  std::vector<reco::Vertex>::const_iterator vertexIt  = vertexes->begin();
  std::vector<reco::Vertex>::const_iterator vertexEnd = vertexes->end();

  for (; vertexIt != vertexEnd; ++vertexIt) 
    {

      const reco::Vertex& vertex = *vertexIt;

      if (!vertex.isValid()) continue;
      ++nVtx;

      if (vertexIt == vertexes->begin()) 
	{
	  event_.primaryVertex[0] = vertex.x();
	  event_.primaryVertex[1] = vertex.y();
	  event_.primaryVertex[2] = vertex.z();

	  for (unsigned int ix=0; ix<3; ++ix) 
	    {
	      for (unsigned int iy=0; iy<3; ++iy) 
		{
		  event_.cov_primaryVertex[ix][iy] = vertex.covariance(ix,iy);
		}
	    }
	}
    }
  
  event_.nVtx = nVtx;
  
}


Int_t MuonPogTreeProducer::fillMuons(const edm::Handle<edm::View<reco::Muon> > & muons,
				     const edm::Handle<std::vector<reco::Vertex> > & vertexes,
				     const edm::Handle<reco::BeamSpot> & beamSpot)
{
  
  edm::View<reco::Muon>::const_iterator muonIt  = muons->begin();
  edm::View<reco::Muon>::const_iterator muonEnd = muons->end();
  
  for (; muonIt != muonEnd; ++muonIt) 
    {
      

      const reco::Muon& mu = (*muonIt);
      
      bool isGlobal      = mu.isGlobalMuon();
      bool isTracker     = mu.isTrackerMuon();
      bool isTrackerArb  = muon::isGoodMuon(mu, muon::TrackerMuonArbitrated); 
      bool isRPC         = mu.isRPCMuon();
      bool isStandAlone  = mu.isStandAloneMuon();
      bool isPF          = mu.isPFMuon();

      bool hasInnerTrack = !mu.innerTrack().isNull();
      bool hasTunePTrack = !mu.tunePMuonBestTrack().isNull();
      bool hasPickyTrack = !mu.pickyTrack().isNull();
      bool hasDytTrack = !mu.dytTrack().isNull();
      bool hasTpfmsTrack = !mu.tpfmsTrack().isNull();
      
      muon_pog::Muon ntupleMu;
      
      ntupleMu.pt     = mu.pt();
      ntupleMu.eta    = mu.eta();
      ntupleMu.phi    = mu.phi();
      ntupleMu.charge = mu.charge();

      ntupleMu.fits.push_back(muon_pog::MuonFit(mu.pt(),mu.eta(),mu.phi(),
						mu.charge(),mu.muonBestTrack()->ptError()));

      ntupleMu.fits.push_back(muon_pog::MuonFit(hasInnerTrack ? mu.innerTrack()->pt()  : -1000.,
      						hasInnerTrack ? mu.innerTrack()->eta() : -1000.,
      						hasInnerTrack ? mu.innerTrack()->phi() : -1000.,
      						hasInnerTrack ? mu.innerTrack()->charge()  : -1000.,
      						hasInnerTrack ? mu.innerTrack()->ptError() : -1000.));

      ntupleMu.fits.push_back(muon_pog::MuonFit(isStandAlone ? mu.outerTrack()->pt()  : -1000.,
      						isStandAlone ? mu.outerTrack()->eta() : -1000.,
      						isStandAlone ? mu.outerTrack()->phi() : -1000.,
      						isStandAlone ? mu.outerTrack()->charge()  : -1000.,
      						isStandAlone ? mu.outerTrack()->ptError() : -1000.));

      ntupleMu.fits.push_back(muon_pog::MuonFit(isGlobal ? mu.globalTrack()->pt()  : -1000.,
      						isGlobal ? mu.globalTrack()->eta() : -1000.,
      						isGlobal ? mu.globalTrack()->phi() : -1000.,
      						isGlobal ? mu.globalTrack()->charge()  : -1000.,
      						isGlobal ? mu.globalTrack()->ptError() : -1000.));

      ntupleMu.fits.push_back(muon_pog::MuonFit(hasTunePTrack ? mu.tunePMuonBestTrack()->pt()  : -1000.,
      						hasTunePTrack ? mu.tunePMuonBestTrack()->eta() : -1000.,
      						hasTunePTrack ? mu.tunePMuonBestTrack()->phi() : -1000.,
      						hasTunePTrack ? mu.tunePMuonBestTrack()->charge()  : -1000.,
      						hasTunePTrack ? mu.tunePMuonBestTrack()->ptError() : -1000.));
      
      ntupleMu.fits.push_back(muon_pog::MuonFit(hasPickyTrack ? mu.pickyTrack()->pt()  : -1000.,
                        hasPickyTrack ? mu.pickyTrack()->eta() : -1000.,
                        hasPickyTrack ? mu.pickyTrack()->phi() : -1000.,
                        hasPickyTrack ? mu.pickyTrack()->charge()  : -1000.,
                        hasPickyTrack ? mu.pickyTrack()->ptError() : -1000.));
      
      ntupleMu.fits.push_back(muon_pog::MuonFit(hasDytTrack ? mu.dytTrack()->pt()  : -1000.,
                        hasDytTrack ? mu.dytTrack()->eta() : -1000.,
                        hasDytTrack ? mu.dytTrack()->phi() : -1000.,
                        hasDytTrack ? mu.dytTrack()->charge()  : -1000.,
                        hasDytTrack ? mu.dytTrack()->ptError() : -1000.));
      
      ntupleMu.fits.push_back(muon_pog::MuonFit(hasTpfmsTrack ? mu.tpfmsTrack()->pt()  : -1000.,
                        hasTpfmsTrack ? mu.tpfmsTrack()->eta() : -1000.,
                        hasTpfmsTrack ? mu.tpfmsTrack()->phi() : -1000.,
                        hasTpfmsTrack ? mu.tpfmsTrack()->charge()  : -1000.,
                        hasTpfmsTrack ? mu.tpfmsTrack()->ptError() : -1000.));

      // Detector Based Isolation
      reco::MuonIsolation detIso03 = mu.isolationR03();

      ntupleMu.trackerIso = detIso03.sumPt;
      ntupleMu.EMCalIso   = detIso03.emEt;
      ntupleMu.HCalIso    = detIso03.hadEt;
      
      // PF Isolation
      reco::MuonPFIsolation pfIso04 = mu.pfIsolationR04();
      reco::MuonPFIsolation pfIso03 = mu.pfIsolationR03();

      ntupleMu.chargedHadronIso   = pfIso04.sumChargedHadronPt;
      ntupleMu.chargedHadronIsoPU = pfIso04.sumPUPt; 
      ntupleMu.neutralHadronIso   = pfIso04.sumNeutralHadronEt;
      ntupleMu.photonIso          = pfIso04.sumPhotonEt;

      ntupleMu.isGlobal     = isGlobal ? 1 : 0;	
      ntupleMu.isTracker    = isTracker ? 1 : 0;	
      ntupleMu.isTrackerArb = isTrackerArb ? 1 : 0;	
      ntupleMu.isRPC        = isRPC ? 1 : 0;
      ntupleMu.isStandAlone = isStandAlone ? 1 : 0;
      ntupleMu.isPF         = isPF ? 1 : 0;

      ntupleMu.algo      = hasInnerTrack ? mu.innerTrack()->algo() : -999;
      ntupleMu.origAlgo  = hasInnerTrack ? mu.innerTrack()->originalAlgo()  : -999;
 
      ntupleMu.nHitsGlobal     = isGlobal     ? mu.globalTrack()->numberOfValidHits() : -999;	
      ntupleMu.nHitsTracker    = isTracker    ? mu.innerTrack()->numberOfValidHits()  : -999;	
      ntupleMu.nHitsStandAlone = isStandAlone ? mu.outerTrack()->numberOfValidHits()  : -999;

      ntupleMu.glbNormChi2              = isGlobal      ? mu.globalTrack()->normalizedChi2() : -999; 
      ntupleMu.trkNormChi2	        = hasInnerTrack ? mu.innerTrack()->normalizedChi2()  : -999; 
      ntupleMu.trkMuonMatchedStations   = isTracker     ? mu.numberOfMatchedStations()       : -999; 
      ntupleMu.glbMuonValidHits	        = isGlobal      ? mu.globalTrack()->hitPattern().numberOfValidMuonHits()       : -999; 
      ntupleMu.trkPixelValidHits	= hasInnerTrack ? mu.innerTrack()->hitPattern().numberOfValidPixelHits()       : -999; 
      ntupleMu.trkPixelLayersWithMeas   = hasInnerTrack ? mu.innerTrack()->hitPattern().pixelLayersWithMeasurement()   : -999; 
      ntupleMu.trkTrackerLayersWithMeas = hasInnerTrack ? mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() : -999; 

      ntupleMu.bestMuPtErr              = mu.muonBestTrack()->ptError(); 

      ntupleMu.trkValidHitFrac = hasInnerTrack           ? mu.innerTrack()->validFraction()       : -999; 
      ntupleMu.trkStaChi2      = isGlobal                ? mu.combinedQuality().chi2LocalPosition : -999; 
      ntupleMu.trkKink         = isGlobal                ? mu.combinedQuality().trkKink           : -999; 
      ntupleMu.muSegmComp      = (isGlobal || isTracker) ? muon::segmentCompatibility(mu)         : -999; 

      ntupleMu.isTrkMuOST               = muon::isGoodMuon(mu, muon::TMOneStationTight) ? 1 : 0; 
      ntupleMu.isTrkHP                  = hasInnerTrack && mu.innerTrack()->quality(reco::TrackBase::highPurity) ? 1 : 0; 

      if ( mu.isMatchesValid() )
        {
          for ( reco::MuonChamberMatch match : mu.matches() )
            {
              muon_pog::ChambMatch ntupleMatch;
	      
              if ( getMuonChamberId(match.id,
                                    ntupleMatch.type,ntupleMatch.id_r,
                                    ntupleMatch.id_phi,ntupleMatch.id_eta)
		   )
                {
		  
                  ntupleMatch.x = mu.trackX(match.station(),match.detector());
                  ntupleMatch.y = mu.trackY(match.station(),match.detector());
		  
                  ntupleMatch.dXdZ = mu.trackDxDz(match.station(),match.detector());
                  ntupleMatch.dYdZ = mu.trackDyDz(match.station(),match.detector());
		  
                  ntupleMatch.errx = mu.trackXErr(match.station(),match.detector());
                  ntupleMatch.erry = mu.trackYErr(match.station(),match.detector());
		  
                  ntupleMatch.errDxDz = mu.trackDxDzErr(match.station(),match.detector());
                  ntupleMatch.errDyDz = mu.trackDyDzErr(match.station(),match.detector());
		  
                  if (ntupleMatch.type == muon_pog::MuonDetType::DT)
                    {
		      
                      bool hasMatch = false;
                      std::size_t iDtSeg = 0;
		      
                      for ( const auto & dtSegment : event_.dtSegments )
                        {
                          for ( const auto & segMatch : match.segmentMatches )
                            {
			      
                              std::size_t refKey = 
				segMatch.dtSegmentRef.isNonnull() ? 
				segMatch.dtSegmentRef.key() : 
				std::numeric_limits<std::size_t>::max();
			      
                              if ( refKey == iDtSeg )
                                {
				  bool isArb = segMatch.isMask(reco::MuonSegmentMatch::BestInChamberByDR) &&
				    segMatch.isMask(reco::MuonSegmentMatch::BelongsToTrackByDR);
				  
				  ntupleMatch.indexes.push_back(refKey);
				  ntupleMatch.matchQuals.push_back(isArb ? muon_pog::MuonMatchQual::TrkArb :
								   muon_pog::MuonMatchQual::Trk); 
				  hasMatch = true;
				  break;
                                }

                            }
			  
			  if ( !hasMatch &&
			       dtSegment.id_r   == ntupleMatch.id_r &&
			       dtSegment.id_phi == ntupleMatch.id_phi &&
			       dtSegment.id_eta == ntupleMatch.id_eta )
                            {
			      ntupleMatch.indexes.push_back(iDtSeg);
			      ntupleMatch.matchQuals.push_back(muon_pog::MuonMatchQual::None);
                            }
                    ++iDtSeg;
                        }
		    }
		
                  if (ntupleMatch.type == muon_pog::MuonDetType::CSC)
                    {
        
                      bool hasMatch = false;
                      std::size_t iCscSeg = 0;      
		      
                      for ( const auto & cscSegment : event_.cscSegments )
                        {
                          for ( const auto & segMatch : match.segmentMatches )
                            {
			      
                              std::size_t refKey = segMatch.cscSegmentRef.isNonnull() ? 
				segMatch.cscSegmentRef.key() : 
				std::numeric_limits<std::size_t>::max();
			      
                              if ( refKey == iCscSeg )
                                {
				  bool isArb = segMatch.isMask(reco::MuonSegmentMatch::BestInChamberByDR) &&
				    segMatch.isMask(reco::MuonSegmentMatch::BelongsToTrackByDR);
				  
				  ntupleMatch.indexes.push_back(refKey);
				  ntupleMatch.matchQuals.push_back(isArb ? muon_pog::MuonMatchQual::TrkArb :
								   muon_pog::MuonMatchQual::Trk); 
				  hasMatch = true;
				  break;
                                }
			      
                            }

			  if ( !hasMatch &&
			       cscSegment.id_r   == ntupleMatch.id_r &&
			       cscSegment.id_phi == ntupleMatch.id_phi &&
			       cscSegment.id_eta == ntupleMatch.id_eta )
                            {
			      ntupleMatch.indexes.push_back(iCscSeg);
			      ntupleMatch.matchQuals.push_back(muon_pog::MuonMatchQual::None);
                            }
                    ++iCscSeg;
                        }
                    }
		  
                }

      ntupleMu.matches.push_back(ntupleMatch); 

            }

	}
          
      
      ntupleMu.dxyBest  = -999; 
      ntupleMu.dzBest   = -999; 
      ntupleMu.dxyInner = -999; 
      ntupleMu.dzInner  = -999; 

      ntupleMu.isoPflow04 = (pfIso04.sumChargedHadronPt+ 
			     std::max(0.,pfIso04.sumPhotonEt+pfIso04.sumNeutralHadronEt - 0.5*pfIso04.sumPUPt)) / mu.pt();
    
      ntupleMu.isoPflow03 = (pfIso03.sumChargedHadronPt+ 
			     std::max(0.,pfIso03.sumPhotonEt+pfIso03.sumNeutralHadronEt - 0.5*pfIso03.sumPUPt)) / mu.pt();

      double dxybs = isGlobal ? mu.globalTrack()->dxy(beamSpot->position()) :
	hasInnerTrack ? mu.innerTrack()->dxy(beamSpot->position()) : -1000;
      double dzbs  = isGlobal ? mu.globalTrack()->dz(beamSpot->position()) :
	hasInnerTrack ? mu.innerTrack()->dz(beamSpot->position()) : -1000;

      double dxy = -1000.;
      double dz  = -1000.;

      ntupleMu.isSoft    = 0;	  
      ntupleMu.isTight   = 0;	  
      ntupleMu.isHighPt  = 0;
      ntupleMu.isLoose   = muon::isLooseMuon(mu)  ? 1 : 0;	  
      ntupleMu.isMedium  = muon::isMediumMuon(mu) ? 1 : 0;	  

      if (vertexes->size() > 0)
	{
	  const reco::Vertex & vertex = vertexes->at(0);

	  dxy = isGlobal ? mu.globalTrack()->dxy(vertex.position()) :
	    hasInnerTrack ? mu.innerTrack()->dxy(vertex.position()) : -1000;
	  dz = isGlobal ? mu.globalTrack()->dz(vertex.position()) :
	    hasInnerTrack ? mu.innerTrack()->dz(vertex.position()) : -1000;
 
	  ntupleMu.dxyBest  = mu.muonBestTrack()->dxy(vertex.position()); 
	  ntupleMu.dzBest   = mu.muonBestTrack()->dz(vertex.position()); 
	  if(hasInnerTrack) { 
	    ntupleMu.dxyInner = mu.innerTrack()->dxy(vertex.position()); 
	    ntupleMu.dzInner  = mu.innerTrack()->dz(vertex.position()); 
	  } 

	  ntupleMu.isSoft    = muon::isSoftMuon(mu,vertex)   ? 1 : 0;	  
	  ntupleMu.isTight   = muon::isTightMuon(mu,vertex)  ? 1 : 0;	  
	  ntupleMu.isHighPt  = muon::isHighPtMuon(mu,vertex) ? 1 : 0;

	}

      ntupleMu.dxy    = dxy;
      ntupleMu.dz     = dz;
      ntupleMu.edxy   = isGlobal ? mu.globalTrack()->dxyError() : hasInnerTrack ? mu.innerTrack()->dxyError() : -1000;
      ntupleMu.edz    = isGlobal ? mu.globalTrack()->dzError()  : hasInnerTrack ? mu.innerTrack()->dzError() : -1000;

      ntupleMu.dxybs  = dxybs;
      ntupleMu.dzbs   = dzbs;

      if(mu.isTimeValid()) { 
	ntupleMu.muonTimeDof = mu.time().nDof; 
	ntupleMu.muonTime    = mu.time().timeAtIpInOut; 
	ntupleMu.muonTimeErr = mu.time().timeAtIpInOutErr; 
      } 
      else { 
	ntupleMu.muonTimeDof = -999; 
	ntupleMu.muonTime    = -999; 
	ntupleMu.muonTimeErr = -999; 
      } 

      if(mu.rpcTime().nDof > 0) { 
	ntupleMu.muonRpcTimeDof = mu.rpcTime().nDof; 
	ntupleMu.muonRpcTime    = mu.rpcTime().timeAtIpInOut; 
	ntupleMu.muonRpcTimeErr = mu.rpcTime().timeAtIpInOutErr; 
      } 
      else { 
	ntupleMu.muonRpcTimeDof = -999; 
	ntupleMu.muonRpcTime    = -999; 
	ntupleMu.muonRpcTimeErr = -999; 
      } 

      // asking for a TRK or GLB muon with minimal pT cut
      // ignoring STA muons in this logic
      if ( m_minMuPtCut < 0 ||
	   (
	    (isTracker || isGlobal || isStandAlone) &&
	    (ntupleMu.fitPt(muon_pog::MuonFitType::DEFAULT) > m_minMuPtCut ||
	     ntupleMu.fitPt(muon_pog::MuonFitType::GLB)     > m_minMuPtCut ||
	     ntupleMu.fitPt(muon_pog::MuonFitType::TUNEP)   > m_minMuPtCut ||
	     ntupleMu.fitPt(muon_pog::MuonFitType::INNER)   > m_minMuPtCut ||
	     ntupleMu.fitPt(muon_pog::MuonFitType::STA)     > m_minMuPtCut ||
	     ntupleMu.fitPt(muon_pog::MuonFitType::PICKY)   > m_minMuPtCut ||
	     ntupleMu.fitPt(muon_pog::MuonFitType::DYT)     > m_minMuPtCut ||
	     ntupleMu.fitPt(muon_pog::MuonFitType::TPFMS)   > m_minMuPtCut)
	     )
          )
      {
        event_.muons.push_back(ntupleMu);
        
        std::vector<Float_t> PTs = {ntupleMu.fitPt(muon_pog::MuonFitType::DEFAULT),
				    ntupleMu.fitPt(muon_pog::MuonFitType::GLB),
				    ntupleMu.fitPt(muon_pog::MuonFitType::TUNEP),
				    ntupleMu.fitPt(muon_pog::MuonFitType::INNER),
				    ntupleMu.fitPt(muon_pog::MuonFitType::PICKY),
				    ntupleMu.fitPt(muon_pog::MuonFitType::DYT),
				    ntupleMu.fitPt(muon_pog::MuonFitType::TPFMS)};
        eventId_.maxPTs.push_back(*std::max_element(PTs.begin(), PTs.end()));
      }

    }

  return event_.muons.size();

}



void MuonPogTreeProducer::fillDtSegments(const edm::Handle<DTRecSegment4DCollection> & segments)
{
  
  DTRecSegment4DCollection::id_iterator chambIt  = segments->id_begin();
  DTRecSegment4DCollection::id_iterator chambEnd = segments->id_end();

  
  for (; chambIt != chambEnd; ++chambIt)
    {
    
      auto range = segments->get(*chambIt);
      for (auto segment = range.first; segment != range.second; ++segment)
        {

	  muon_pog::MuonSegment ntupleSeg;
	  muon_pog::MuonDetType dummyType = muon_pog::MuonDetType::RPC;

	  getMuonChamberId((*chambIt),
			   dummyType,ntupleSeg.id_r,
			   ntupleSeg.id_phi,ntupleSeg.id_eta);

	  std::cout << "[MuonPogTreeProducer] : DT segment type : " << dummyType << std::endl;
	  
	  ntupleSeg.x = segment->localPosition().x(); 
	  ntupleSeg.y = segment->localPosition().y();
	  
	  ntupleSeg.phi = 0.; 
	  ntupleSeg.eta = 0.;

	  ntupleSeg.dXdZ = segment->localDirection().x(); 
	  ntupleSeg.dYdZ = segment->localDirection().y();

	  ntupleSeg.errx = 0.; 
	  ntupleSeg.erry = 0.; 

	  ntupleSeg.errDxDz = 0.; 
	  ntupleSeg.errDyDz = 0.; 

	  ntupleSeg.chi2 = segment->hasPhi() ? (segment->phiSegment()->chi2() /
                                            segment->phiSegment()->degreesOfFreedom()) : -999.;  
	  
	  ntupleSeg.time = segment->hasPhi() ? segment->phiSegment()->t0() : -999.; 
	  
	  ntupleSeg.nHitsX = segment->hasPhi() ? segment->phiSegment()->specificRecHits().size() : 0; 
	  ntupleSeg.nHitsY = segment->hasZed() ? segment->zSegment()->specificRecHits().size()   : 0;
	  
	  event_.dtSegments.push_back(ntupleSeg);
	  
        }
    }

}



void MuonPogTreeProducer::fillCscSegments(const edm::Handle<CSCSegmentCollection> & segments)
{
  
  CSCSegmentCollection::id_iterator chambIt  = segments->id_begin();
  CSCSegmentCollection::id_iterator chambEnd = segments->id_end();

  for (; chambIt != chambEnd; ++chambIt)
    {
    
      auto range = segments->get(*chambIt);
      for (auto segment = range.first; segment != range.second; ++segment)
        {

	  muon_pog::MuonSegment ntupleSeg;
	  muon_pog::MuonDetType dummyType = muon_pog::MuonDetType::RPC;
	  
	  getMuonChamberId((*chambIt),
			   dummyType,ntupleSeg.id_r,
			   ntupleSeg.id_phi,ntupleSeg.id_eta);
	  
	  std::cout << "[MuonPogTreeProducer] : CSC segment type : " << dummyType << std::endl;

	  ntupleSeg.x = segment->localPosition().x(); 
	  ntupleSeg.y = segment->localPosition().y();

	  ntupleSeg.phi = 0.; 
	  ntupleSeg.eta = 0.;

	  ntupleSeg.dXdZ = segment->localDirection().x(); 
	  ntupleSeg.dYdZ = segment->localDirection().y();

	  ntupleSeg.errx = 0.; 
	  ntupleSeg.erry = 0.; 

	  ntupleSeg.errDxDz = 0.; 
	  ntupleSeg.errDyDz = 0.; 

	  ntupleSeg.chi2 = segment->chi2() /
                       segment->degreesOfFreedom();  

	  ntupleSeg.time = segment->time(); 
	  
	  ntupleSeg.nHitsX = segment->nRecHits(); 
	  ntupleSeg.nHitsY = 0;

	  event_.cscSegments.push_back(ntupleSeg);
	  
        }
    }
  
}


bool MuonPogTreeProducer::getMuonChamberId(const DetId & id, muon_pog::MuonDetType & det,
					   Int_t & r, Int_t & phi, Int_t & eta) const
{

  if (id.det() == DetId::Muon && id.subdetId() == MuonSubdetId::DT)
    {
      DTChamberId dtId(id.rawId());  
  
      det = muon_pog::MuonDetType::DT;
      r   = dtId.station();
      phi = dtId.sector();
      eta = dtId.wheel();

      return true;
    }

  if (id.det() == DetId::Muon && id.subdetId() == MuonSubdetId::CSC)
    {
      CSCDetId cscId(id.rawId());
    
      det = muon_pog::MuonDetType::CSC;
      r   = cscId.station() * cscId.zendcap();
      phi = cscId.chamber();
      eta = cscId.ring();

      return true;
    }

  return false;
      
}


#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(MuonPogTreeProducer);
