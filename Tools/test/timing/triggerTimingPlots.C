#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "THStack.h"
#include "TLegend.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLorentzVector.h"

#include "interface/MuonPogTree.h"
#include "interface/Utils.h"
#include "tdrstyle.C"

#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <fstream> 
#include <vector>
#include <regex>
#include <map>

#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/exceptions.hpp>
#include <boost/lexical_cast.hpp>

// Helper classes defintion *****
// 1. SampleConfig : configuration class containing sample information
// 2. TagAndProbeConfig : configuration class containing TnP cuts information
// 3. Plotter : class containing the plot definition and defining the plot filling 
//              for a given sample <= CB modify this to add new variables
// ******************************

namespace muon_pog {
 
  class SampleConfig {

  public :

    // config parameters (public for direct access)

    TString fileName;  
    TString sampleName;  
    Float_t cSection;
    Float_t nEvents;
    Bool_t applyReweighting;
    std::vector<int> runs;
        
    SampleConfig() {};
    
#ifndef __MAKECINT__ // CB CINT doesn't like boost :'-(    
    SampleConfig(boost::property_tree::ptree::value_type & vt); 
#endif

    ~SampleConfig() {};

  private:
    std::vector<int> toArray(const std::string & entries); 
    
  };

  class TagAndProbeConfig {

  public :
    
    // config parameters (public for direct access)
    
    Float_t pair_minInvMass;
    Float_t pair_maxInvMass;
    Float_t pair_minCosAngle;      
    Float_t pair_maxPtBalance;      
    Float_t pair_maxDz;      
    Float_t pair_minDr;      

    Float_t tag_minPt;      

    TString     tag_ID;
    Float_t     tag_isoCut;
    Float_t     tag_isoCutAbs;
    Float_t     tag_hltDrCut;
    std::string tag_hltFilter;
    
    std::string muon_trackType; // applies to both, tag and probe     
  
    Float_t probe_minPt;      
    Float_t probe_isoCut;
    Float_t probe_isoCutAbs;
    Float_t probe_maxPrimDphi;      
    Int_t   probe_minPrimBX;      

    Float_t probe_minEtaBX;      
    Float_t probe_maxEtaBX;      


    std::vector<TString> probe_IDs;
    std::vector<std::pair<TString,TString> > probe_fEtaBins;
     
    std::string hlt_path; 
   
    TagAndProbeConfig() {};
    
#ifndef __MAKECINT__ // CB CINT doesn't like boost :'-(    
    TagAndProbeConfig(boost::property_tree::ptree::value_type & vt); 
#endif

    ~TagAndProbeConfig() {};
    
  private:
    std::vector<TString> toArray(const std::string & entries); 
    std::vector<std::pair<TString,TString> > toPairArray(const std::vector<TString> &,
							 const std::vector<TString> &); 
  
  };

  class Plotter {

  public :

    enum HistoType { KIN=0, CONT, TRIG, TIMING};
      
    Plotter(muon_pog::TagAndProbeConfig tnpConfig, muon_pog::SampleConfig & sampleConfig) :
      m_tnpConfig(tnpConfig) , m_sampleConfig(sampleConfig) {};
    ~Plotter() {};
    
    void book(TFile *outFile);
    void fill(const std::vector<muon_pog::Muon> & muons, const muon_pog::HLT & hlt, const muon_pog::Event & ev, float weight);

    std::map<Plotter::HistoType, std::map<TString, TH1* > > m_histos;
    std::map<Plotter::HistoType, std::map<TString, TEfficiency *> > m_effs;
    TagAndProbeConfig m_tnpConfig;
    SampleConfig m_sampleConfig;
        
  };
  
}

// Helper classes defintion *****
// 1. parseConfig : parse the full cfg file
// 1. comparisonPlot : make a plot overlayng data and MC for a given plot
// ******************************

namespace muon_pog {
  void parseConfig(const std::string configFile, TagAndProbeConfig & tpConfig,
		   std::vector<SampleConfig> & sampleConfigs);
  
  void addUnderFlow(TH1 &hist);
  void addOverFlow(TH1 &hist);
}



// The main program******** *****
// 1. Get configuration file and produces configurations
// 2. Create Plotters and loop on the event to fill them
// 3. Writes results in cnfigurable outuput file
// ******************************

int main(int argc, char* argv[]){
  using namespace muon_pog;


  if (argc != 3) 
    {
      std::cout << "Usage : "
		<< argv[0] << " PATH_TO_CONFIG_FILE PATH_TO_OUTPUT_DIR\n";
      exit(100);
    }

  std::string configFile(argv[1]);
  
  std::cout << "[" << argv[0] << "] Using config file " << configFile << std::endl;

  // Output directory
  TString dirName = argv[2];
  system("mkdir -p " + dirName);
  TFile* outputFile = TFile::Open(dirName + "/results.root","RECREATE"); // CB find a better name for output file  

  // Set it to kTRUE if you do not run interactively
  gROOT->SetBatch(kTRUE); 

  // Initialize Root application
  TRint* app = new TRint("CMS Root Application", &argc, argv);

  setTDRStyle();
 
  TagAndProbeConfig tnpConfig;
  std::vector<SampleConfig> sampleConfigs;
  
  parseConfig(configFile,tnpConfig,sampleConfigs);

  std::vector<Plotter> plotters;

  for (auto sampleConfig : sampleConfigs)
    {

      Plotter plotter(tnpConfig, sampleConfig);
      plotter.book(outputFile);
      
      plotters.push_back(plotter);
    }
 
  for (auto plotter : plotters)
    {

      TString fileName = plotter.m_sampleConfig.fileName;
      std::cout << "[" << argv[0] << "] Processing file "
		<< fileName.Data() << std::endl;  
  
      // Initialize pointers to summary and full event structure

      muon_pog::Event*   ev   = new muon_pog::Event();

      TTree* tree;
      TBranch* evBranch;

      // Open file, get tree, set branches

      TFile* inputFile = TFile::Open(fileName,"READONLY");
      tree = (TTree*)inputFile->Get("MUONPOGTREE");
      if (!tree) inputFile->GetObject("MuonPogTree/MUONPOGTREE",tree);

      evBranch = tree->GetBranch("event");
      evBranch->SetAddress(&ev);

      // Watch number of entries
      int nEntries = plotter.m_sampleConfig.nEvents > 0 ? plotter.m_sampleConfig.nEvents : tree->GetEntriesFast();
      std::cout << "[" << argv[0] << "] Number of entries = " << nEntries << std::endl;

      int nFilteredEvents = 0;

      for (Long64_t iEvent=0; iEvent<nEntries; ++iEvent) 
	{
	  if (tree->LoadTree(iEvent)<0) break;

	  if (iEvent % 25000 == 0 )
	    std::cout << "[" << argv[0] << "] processing event : " << iEvent << "\r" << std::flush;
	    
	  evBranch->GetEntry(iEvent);
	  float weight = ev->genInfos.size() > 0 ?
	    ev->genInfos[0].genWeight/fabs(ev->genInfos[0].genWeight) : 1.;
	 
	  // CB to be fixed when kown what to do for MC !!!

	  //	  if(plotter.m_sampleConfig.applyReweighting==true)
	  //  weight *= ev->nVtx < 60 ? PUweight[ev->nVtx] : 0;
	  
	  plotter.fill(ev->muons, ev->hlt, (*ev), weight);
	}
      
      delete ev;
      inputFile->Close();
      std::cout << std::endl;
	   
    }
  
  outputFile->Write();
  
  if (!gROOT->IsBatch()) app->Run();

  return 0;

}


muon_pog::TagAndProbeConfig::TagAndProbeConfig(boost::property_tree::ptree::value_type & vt)
{

  try
    {

      hlt_path = vt.second.get<std::string>("hlt_path");

      pair_minInvMass  = vt.second.get<Float_t>("pair_minInvMass");
      pair_maxInvMass  = vt.second.get<Float_t>("pair_maxInvMass");
      
      pair_minCosAngle  = vt.second.get<Float_t>("pair_minCosAngle");
      pair_maxPtBalance = vt.second.get<Float_t>("pair_maxPtBalance");
      pair_minDr  = vt.second.get<Float_t>("pair_minDr");
      pair_maxDz  = vt.second.get<Float_t>("pair_maxDz");

      tag_minPt     = vt.second.get<Float_t>("tag_minPt");
      tag_ID        = vt.second.get<std::string>("tag_muonID");
      tag_isoCutAbs = vt.second.get<Float_t>("tag_isoCutAbs");
      tag_isoCut    = vt.second.get<Float_t>("tag_isoCut");

      tag_hltFilter = vt.second.get<std::string>("tag_hltFilter");
      tag_hltDrCut  = vt.second.get<Float_t>("tag_hltDrCut");
      
      muon_trackType = vt.second.get<std::string>("muon_trackType");

      probe_minPt     = vt.second.get<Float_t>("probe_minPt");
      probe_isoCut    = vt.second.get<Float_t>("probe_isoCut");
      probe_isoCutAbs = vt.second.get<Float_t>("probe_isoCutAbs");

      probe_minEtaBX  = vt.second.get<Float_t>("probe_minEtaBX");
      probe_maxEtaBX  = vt.second.get<Float_t>("probe_maxEtaBX");
        
      probe_maxPrimDphi = vt.second.get<Float_t>("probe_maxPrimDphi");
      probe_minPrimBX   = vt.second.get<Int_t>("probe_minPrimBX");

      probe_IDs    = toArray(vt.second.get<std::string>("probe_muonIDs"));
      probe_fEtaBins = toPairArray(toArray(vt.second.get<std::string>("probe_fEtaMin")),
				   toArray(vt.second.get<std::string>("probe_fEtaMax")));
  
    }

  catch (boost::property_tree::ptree_bad_data bd)
    {
      std::cout << "[TagAndProbeConfig] Can't get data : has error : "
		<< bd.what() << std::endl;
      throw std::runtime_error("Bad INI variables");
    }

}

muon_pog::SampleConfig::SampleConfig(boost::property_tree::ptree::value_type & vt)
{
 
 try
    {
      fileName     = TString(vt.second.get<std::string>("fileName").c_str());
      sampleName   = TString(vt.first.c_str());
      cSection = vt.second.get<Float_t>("cSection");
      nEvents = vt.second.get<Float_t>("nEvents");
      applyReweighting = vt.second.get<Bool_t>("applyReweighting");
      runs = toArray(vt.second.get<std::string>("runs"));
    }
  
  catch (boost::property_tree::ptree_bad_data bd)
    {
      std::cout << "[TagAndProbeConfig] Can't get data : has error : "
		<< bd.what() << std::endl;
      throw std::runtime_error("Bad INI variables");
    }
  
}

std::vector<int> muon_pog::SampleConfig::toArray(const std::string& entries)
{
  
  std::vector<int> result;
  std::stringstream sentries(entries);
  std::string item;
  while(std::getline(sentries, item, ','))
    result.push_back(atoi(item.c_str()));
  return result;

}


std::vector<TString> muon_pog::TagAndProbeConfig::toArray(const std::string& entries)
{
  
  std::vector<TString> result;
  std::stringstream sentries(entries);
  std::string item;
  while(std::getline(sentries, item, ','))
    result.push_back(TString(item));
  return result;

}

std::vector<std::pair<TString,TString> > muon_pog::TagAndProbeConfig::toPairArray(const std::vector<TString> & fEtaMin,
										  const std::vector<TString> & fEtaMax)
{

  std::vector<std::pair<TString,TString> > result;

  std::vector<TString>::const_iterator fEtaMinIt  = fEtaMin.begin();
  std::vector<TString>::const_iterator fEtaMinEnd = fEtaMin.end();

  std::vector<TString>::const_iterator fEtaMaxIt  = fEtaMax.begin();
  std::vector<TString>::const_iterator fEtaMaxEnd = fEtaMax.end();
  
  for (; fEtaMinIt != fEtaMinEnd || fEtaMaxIt != fEtaMaxEnd; ++fEtaMinIt, ++fEtaMaxIt)
    {
      result.push_back(std::make_pair((*fEtaMinIt),
				      (*fEtaMaxIt)));
    }

  return result;

}


void muon_pog::Plotter::book(TFile *outFile)
{

  TString sampleTag = m_sampleConfig.sampleName;
  
  outFile->cd("/");
  outFile->mkdir(sampleTag);
  outFile->cd(sampleTag);

  outFile->mkdir(sampleTag+"/efficiencies");      
  outFile->mkdir(sampleTag+"/kinematical_variables");
  outFile->mkdir(sampleTag+"/control");
  outFile->mkdir(sampleTag+"/trigger");
  outFile->mkdir(sampleTag+"/timing");
  
  for (auto fEtaBin : m_tnpConfig.probe_fEtaBins)
    {
      TString etaTag = "_fEtaMin" + fEtaBin.first + "_fEtaMax" + fEtaBin.second;
      
      outFile->cd(sampleTag+"/control");
      
      m_histos[CONT]["invMass" + etaTag] = new TH1F("invMass" + etaTag + sampleTag , 
						    "invMass" + etaTag + sampleTag+ ";mass (GeV); # entries", 
						    200, 0.,1000.);
      m_histos[CONT]["dZ" + etaTag] = new TH1F("dZ" + etaTag + sampleTag , 
					       "dZ" + etaTag + sampleTag+ ";dZ(tag,probe); # entries", 
						    200, -1.,1.);
      m_histos[CONT]["tagPtVsProbePt" + etaTag] = new TH2F("tagPtVsProbePt_" + etaTag, 
							   "tagPtVsProbePt_" + etaTag +";p_{T} tag (GeV/c);p_{T} probe (GeV/c)", 
							   300, 0., 1500., 300, 0., 1500.);
      m_histos[CONT]["probePtVsDr" + etaTag] = new TH2F("probePtVsDr_" + etaTag, 
							"probePtVsDr_" + etaTag +";p_{T} probe (GeV/c);dR(tag,probe)", 
							300, 0., 1500., 72, 0., 2*TMath::Pi());
      m_histos[CONT]["tagEtaVsProbeEta" + etaTag] = new TH2F("tagEtaVsProbeEta_" + etaTag, 
							     "tagEtaVsProbeEta_" + etaTag +";#eta_{tag};#eta_{probe}", 
							     48., -2.4, 2.4, 24., -1.2, 1.2);

      for ( auto & probe_ID : m_tnpConfig.probe_IDs)
	{
	  TString IDTag = "_" + probe_ID;
	  
	  const Double_t ptBins[18] = {0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.,
	                               120., 140, 200., 400., 600., 800., 1200.};

	  TString completeTag = etaTag + IDTag + sampleTag;	  
	  
	  for (Int_t iChamb = 1; iChamb<=4; ++iChamb)
	    {

	      TString chTag(std::to_string(iChamb).c_str());	  

	      outFile->cd(sampleTag+"/control");

	      m_histos[CONT]["sectorVsWheelMB" + chTag + etaTag + IDTag] = new TH2F("sectorVsWheelMB" + chTag + completeTag, 
										    "sectorVsWheelMB" + chTag + completeTag +
										    "sector;wheel", 
										    14, 0.5, 14.5, 5, -2.5, 2.5);

	      outFile->cd(sampleTag+"/trigger");

	      m_histos[TRIG]["nTrigMB" + chTag + etaTag + IDTag]     = new TH1F("nTrigMB" + chTag + completeTag,
										"nTrigMB" + chTag + completeTag +
										";# matched primitives ;# entries",
										15,-.5,14.5);

	      m_histos[TRIG]["qualTrigMB" + chTag + etaTag + IDTag]  = new TH1F("qualTrigMB" + chTag + completeTag,
										"qualTrigMB" + chTag + completeTag +
										";quality;# entries",
										7,-.5,6.5);
	      
	      m_histos[TRIG]["dPhiMB" + chTag + etaTag + IDTag]      = new TH1F("dPhiMB" + chTag + completeTag,
										"dPhiMB" + chTag + completeTag +
										";dPhi(track,primitive);# entries",
										300,-0.5,0.5);

	      m_histos[TRIG]["dPhiVsEtaMB" + chTag + etaTag + IDTag] = new TH2F("dPhiVsEtaMB" + chTag + completeTag,
										"dPhiVsEtaMB" + chTag + completeTag +
										";dPhi(track,primitive);#eta",
										300,-0.5,0.5, 24., -1.2, 1.2);

	      m_histos[TRIG]["dPhiVsPhiMB" + chTag + etaTag + IDTag] = new TH2F("dPhiVsPhiMB" + chTag + completeTag,
										"dPhiVsPhiMB" + chTag + completeTag +
										";dPhi(track,primitive);#phi",
										300,-0.5,0.5, 48,-TMath::Pi(),TMath::Pi());

	      m_histos[TRIG]["dPhiVsPtMB" + chTag + etaTag + IDTag]  = new TH2F("dPhiVsPtMB" + chTag + completeTag,
										"dPhivsPtMB" + chTag + completeTag +
										";dPhi(track,primitive);p_{T}",
										300,-0.5,0.5, 17,ptBins);

	      outFile->cd(sampleTag+"/efficiencies");

	      m_effs[TIMING]["earlyEffVsPtMB" + chTag + etaTag + IDTag] = new TEfficiency("earlyEffVsPtMB" + chTag + completeTag,
											  "earlyEffVsPtMB" + chTag + completeTag +
											  ";p_{T} (GeV/c); fraction of prefiring",
											  17,ptBins);

	      m_effs[TIMING]["earlyEffVsEtaMB" + chTag + "0" + etaTag + IDTag] = new TEfficiency("earlyEffVsEtaMB" + chTag + "0" + completeTag,
												 "earlyEffVsEtaMB" + chTag + "0" + completeTag +
												 ";p_{T} (GeV/c); fraction of prefiring",
												 24., -0.9, 0.9);

	      m_effs[TIMING]["earlyEffVsEtaMB" + chTag + "1" + etaTag + IDTag] = new TEfficiency("earlyEffVsEtaMB" + chTag + "1" + completeTag,
												 "earlyEffVsEtaMB" + chTag + "1" + completeTag +
												 ";p_{T} (GeV/c); fraction of prefiring",
												 24., -0.9, 0.9);

	      m_effs[TIMING]["lateEffVsPtMB" + chTag + etaTag + IDTag] = new TEfficiency("lateEffVsPtMB" + chTag + completeTag,
											 "lateEffVsPtMB" + chTag + completeTag +
											 ";p_{T} (GeV/c); fraction of prefiring",
											 17,ptBins);

	      m_effs[TIMING]["lateEffVsEtaMB" + chTag + "0" + etaTag + IDTag] = new TEfficiency("lateEffVsEtaMB" + chTag + "0" + completeTag,
												"lateEffVsEtaMB" + chTag + "0" + completeTag +
												";p_{T} (GeV/c); fraction of prefiring",
												24., -0.9, 0.9);

	      m_effs[TIMING]["lateEffVsEtaMB" + chTag + "1" + etaTag + IDTag] = new TEfficiency("lateEffVsEtaMB" + chTag + "1" + completeTag,
												"lateEffVsEtaMB" + chTag + "1" + completeTag +
												";p_{T} (GeV/c); fraction of prefiring",
												24., -0.9, 0.9);
	      
	      m_effs[TIMING]["bxm1EffVsPtMB" + chTag + etaTag + IDTag] = new TEfficiency("bxm1EffVsPtMB" + chTag + completeTag,
											 "bxm1EffVsPtMB" + chTag + completeTag +
											 ";p_{T} (GeV/c); fraction of prefiring (BX = -1)",
											 17,ptBins);

	      m_effs[TIMING]["bxm2EffVsPtMB" + chTag + etaTag + IDTag] = new TEfficiency("bxm2EffVsPtMB" + chTag + completeTag,
											 "bxm2EffVsPtMB" + chTag + completeTag +
											 ";p_{T} (GeV/c); fraction of prefiring (BX = -2)",
											 17,ptBins);
	      
	      m_effs[TIMING]["earlyEffVsPhiMB" + chTag + etaTag + IDTag] = new TEfficiency("earlyEffVsPhiMB" + chTag + completeTag,
											   "earlyEffVsPhiMB" + chTag + completeTag +
											   ";p_{T} (GeV/c); fraction of prefiring",
											   48,-TMath::Pi(),TMath::Pi());

	      m_effs[TIMING]["lateEffVsPhiMB" + chTag + etaTag + IDTag] = new TEfficiency("lateEffVsPhiMB" + chTag + completeTag,
											  "lateEffVsPhiMB" + chTag + completeTag +
											  ";p_{T} (GeV/c); fraction of prefiring",
											  48,-TMath::Pi(),TMath::Pi());

	      m_effs[TIMING]["firstEffVsPtMB" + chTag + etaTag + IDTag]  = new TEfficiency("firstEffVsPtMB" + chTag + completeTag,
											   "firstEffVsPtMB" + chTag + completeTag +
											   ";p_{T} (GeV/c); fraction of prefiring",
											   17,ptBins);

	      outFile->cd(sampleTag+"/timing");

	      m_histos[TIMING]["bxTrigMB" + chTag + etaTag + IDTag] = new TH1F("bxTrigMB" + chTag + completeTag,
									       "bxTrigMB" + chTag + completeTag +
									       ";bxTrig(track,primitive);# entries",
									       5,-2.5,2.5);

	      m_histos[TIMING]["bxTrigVsPtMB" + chTag + etaTag + IDTag] = new TProfile("bxTrigVsPtMB" + chTag + completeTag,
										       "bxTrigVsPtMB" + chTag + completeTag +
										       ";bxTrig(track,primitive);#eta",
										       17,ptBins, -2.5, 2.5);

	      m_histos[TIMING]["highPtTrigVsPhiMB" + chTag + etaTag + IDTag] = new TH1F("highPtTrigVsPhiMB" + chTag + completeTag,
											"highPtTrigVsPhiMB" + chTag + completeTag +
											";earlyTrig(track,primitive);#phi",
											48, -TMath::Pi(),TMath::Pi());

	      m_histos[TIMING]["highPtTrigBxMB" + chTag + etaTag + IDTag]  = new TH1F("highPtTrigBxMB" + chTag + completeTag,
										      "highPtTrigBxMB" + chTag + completeTag +
										      ";BX;# entries",
										      5, -2.5, 2.5);

	      m_histos[TIMING]["bxTrigVsEtaMB" + chTag + "0" + etaTag + IDTag]  = new TProfile("bxTrigVsEtaMB" + chTag + "0" + completeTag,
												"bxTrigVsEtaMB" + chTag + "0" + completeTag +
												";#eta;<BX>",
												36, -0.9, 0.9, -2.5, 2.5);

	      m_histos[TIMING]["bxTrigVsEtaMB" + chTag + "1" + etaTag + IDTag]  = new TProfile("bxTrigVsEtaMB" + chTag + "1" + completeTag,
												"bxTrigVsEtaMB" + chTag + "1" + completeTag +
												";#eta;<BX>",
												36, -0.9, 0.9, -2.5, 2.5);
	    }

	  outFile->cd(sampleTag+"/kinematical_variables");

	  m_histos[KIN]["ProbePt" + etaTag + IDTag]  = new TH1F("hProbePt" + completeTag,
								"hProbePt" + completeTag +
								";p_{T} (GeV);# entries",
								300,0.,1500.);
	  m_histos[KIN]["ProbeEta" + etaTag + IDTag] = new TH1F("hProbeEta" + completeTag,
								"hProbeEta" + completeTag +
								";#eta;# entries",
								48,-2.4, 2.4);
	  m_histos[KIN]["ProbePhi" + etaTag + IDTag] = new TH1F("hProbePhi" + completeTag,
								"hProbePhi" + completeTag +
								";#phi;# entries",
								48,-TMath::Pi(),TMath::Pi()); 
	  
	}

    }
  
}

void muon_pog::Plotter::fill(const std::vector<muon_pog::Muon> & muons,
			     const muon_pog::HLT & hlt, const muon_pog::Event & ev, float weight)
{

  bool isGoodRun = false;

  for (auto run :m_sampleConfig.runs)
    {
      if (run == 0 || run == ev.runNumber)
	{
	  isGoodRun = true;
	  break;
	}
    }

  if (!isGoodRun) return;
  
  if (!muon_pog::pathHasFired(hlt,m_tnpConfig.hlt_path)) return;

  std::vector<const muon_pog::Muon *> tagMuons;

  for (auto & tagMuon : muons)
    {

      TLorentzVector tagMuTk(muon_pog::muonTk(tagMuon,m_tnpConfig.muon_trackType));

      if ( tagMuTk.Pt() > m_tnpConfig.tag_minPt                       &&
	   muon_pog::hasFilterMatch(tagMuon,hlt,
				    m_tnpConfig.tag_hltFilter,
				    m_tnpConfig.tag_hltDrCut)         &&
	   muon_pog::hasGoodId(tagMuon,m_tnpConfig.tag_ID)            && 
	   tagMuon.trackerIso / tagMuTk.Pt() < m_tnpConfig.tag_isoCut &&
	   tagMuon.trackerIso < m_tnpConfig.tag_isoCutAbs
	 )
	tagMuons.push_back(&tagMuon);
    }
  
  std::vector<const muon_pog::Muon *> probeMuons;
  
  for (auto tagMuonPointer : tagMuons)
    {
      const muon_pog::Muon & tagMuon = *tagMuonPointer;

      bool hasGoodPair = false;
      
      for (auto & muon : muons)
	{
	  
	  if ( tagMuonPointer != &muon && 
	       muon_pog::chargeFromTrk(tagMuon,m_tnpConfig.muon_trackType) *
	       muon_pog::chargeFromTrk(muon,m_tnpConfig.muon_trackType) == -1)
	    {

	      TLorentzVector tagMuTk(muon_pog::muonTk(tagMuon,m_tnpConfig.muon_trackType));
	      TLorentzVector muTk(muon_pog::muonTk(muon,m_tnpConfig.muon_trackType));

	      bool hasGoodId = false;
	      
	      for (auto & probe_ID : m_tnpConfig.probe_IDs)
		{		  
		  if(muon_pog::hasGoodId(muon,probe_ID)) 
		    {
		      hasGoodId = true;
		      break;		      
		    }
		}

	      // General Probe Muons	      
	      if(hasGoodId                            && 
		 muTk.Pt() > m_tnpConfig.probe_minPt  &&
		 muon.trackerIso                 < m_tnpConfig.probe_isoCutAbs   &&
		 muon.trackerIso / muTk.Pt()     < m_tnpConfig.probe_isoCut      &&
		 pTBalance(tagMuTk,muTk)         < m_tnpConfig.pair_maxPtBalance &&
		 pairCos3DAngle(tagMuTk,muTk)    > m_tnpConfig.pair_minCosAngle  &&
		 tagMuTk.DeltaR(muTk)            > m_tnpConfig.pair_minDr   &&
		 std::abs(pairDZ(tagMuon,muon))  < m_tnpConfig.pair_maxDz)
		{
		  
		  Float_t mass = (tagMuTk+muTk).M();
		  Float_t dilepPt = (tagMuTk+muTk).Pt();
		  
		  for (auto fEtaBin : m_tnpConfig.probe_fEtaBins)
		    {
		      if (fabs(muTk.Eta()) > fEtaBin.first.Atof() &&
			  fabs(muTk.Eta()) < fEtaBin.second.Atof() )
			{
			  
			  TString etaTag = "_fEtaMin" + fEtaBin.first + "_fEtaMax" + fEtaBin.second;		  	  

			  m_histos[CONT]["invMass" +  etaTag]->Fill(mass,weight);

			  if ( mass > m_tnpConfig.pair_minInvMass & 
			       mass < m_tnpConfig.pair_maxInvMass )
			    {

			      m_histos[CONT]["dZ" + etaTag]->Fill(pairDZ(tagMuon,muon), weight);
			      static_cast<TH2F*>(m_histos[CONT]["tagPtVsProbePt" + etaTag])->Fill(tagMuTk.Pt(),muTk.Pt(), weight);
			      static_cast<TH2F*>(m_histos[CONT]["probePtVsDr" + etaTag])->Fill(muTk.Pt(), tagMuTk.DeltaR(muTk), weight);
			      static_cast<TH2F*>(m_histos[CONT]["tagEtaVsProbeEta" + etaTag])->Fill(tagMuTk.Eta(),muTk.Eta(), weight);
			      probeMuons.push_back(&muon);
			      hasGoodPair = true;
			      break; // CB If a muon is already a probe don't loop on other tags
			      
			    }
			}
		    }
		}

	    }

	  if(hasGoodPair)
	    break;

	}

      
    }
  
  for (auto probeMuonPointer : probeMuons)
    {
      const muon_pog::Muon & probeMuon = *probeMuonPointer;
      
      TLorentzVector probeMuTk(muonTk(probeMuon,m_tnpConfig.muon_trackType));

      for (auto fEtaBin : m_tnpConfig.probe_fEtaBins)
	{
	  
	  if (fabs(probeMuTk.Eta()) > fEtaBin.first.Atof() &&
	      fabs(probeMuTk.Eta()) < fEtaBin.second.Atof() )
	    {
	      
	      TString etaTag = "_fEtaMin" + fEtaBin.first + "_fEtaMax" + fEtaBin.second;
	      
	      for (auto & probe_ID : m_tnpConfig.probe_IDs)
		{
		  TString IDTag = "_" + probe_ID;
		  
		  if(muon_pog::hasGoodId(probeMuon,probe_ID)) 
		    {	 
		      m_histos[KIN]["ProbePt" + etaTag + IDTag]->Fill(probeMuTk.Pt(), weight);
		      m_histos[KIN]["ProbeEta" + etaTag + IDTag]->Fill(probeMuTk.Eta(), weight);
		      m_histos[KIN]["ProbePhi" + etaTag + IDTag]->Fill(probeMuTk.Phi(), weight);
		      
		      Int_t nTrig[4]    = { 0, 0, 0, 0 };
		      Int_t firstBX[4]  = { 9, 9, 9, 9 };
		      
		      for(auto match : probeMuon.matches) 
			{

			  // CB protetction to avoid standalone induced matchings
			  if (match.errx < 0)
			    continue;
			  
			  TString chTag(std::to_string(match.id_r));
			  
			  for (auto index : match.trigIndexes)
			    {

			      Int_t wh0Tag = (match.id_phi%4) < 2 ? 0 : 1;
			      
 			      const auto & dtPrim = ev.dtPrimitives.at(index);
			      Int_t bx = dtPrim.bxTrackFinder();
			      
			      if (match.id_r   != dtPrim.id_r ||
				  match.id_eta != dtPrim.id_eta ||
				  !( match.id_phi == dtPrim.id_phi ||
				     ( match.id_phi == 13 && dtPrim.id_phi == 4 ) ||
				     ( match.id_phi == 14 && dtPrim.id_phi == 10)
				   ) ||
				  match.type != muon_pog::MuonDetType::DT)
				{
				  std::cout << "CAZZO" << std::endl;
				  continue;
				}

			      Float_t dPhi = match.phi - dtPrim.phiGlb();
			      dPhi = dPhi >  TMath::Pi() ? dPhi - 2*TMath::Pi() :
				     dPhi < -TMath::Pi() ? dPhi + 2*TMath::Pi() :  dPhi;
			    
			      m_histos[TRIG]["dPhiMB" + chTag + etaTag + IDTag]->Fill(dPhi, weight);
			      static_cast<TH2F*>(m_histos[TRIG]["dPhiVsEtaMB" + chTag + etaTag + IDTag])->Fill(dPhi, probeMuTk.Eta(), weight);
			      static_cast<TH2F*>(m_histos[TRIG]["dPhiVsPhiMB" + chTag + etaTag + IDTag])->Fill(dPhi, probeMuTk.Phi(), weight);
			      static_cast<TH2F*>(m_histos[TRIG]["dPhiVsPtMB" + chTag + etaTag + IDTag])->Fill(dPhi, probeMuTk.Pt(),  weight);


			      if (std::abs(dPhi) < m_tnpConfig.probe_maxPrimDphi)
				{ 

				  nTrig[match.id_r - 1]++;
				  static_cast<TH2F*>(m_histos[CONT]["sectorVsWheelMB" + chTag + etaTag + IDTag])->Fill(match.id_phi, match.id_eta, weight);
				  
				  static_cast<TProfile*>(m_histos[TIMING]["bxTrigVsEtaMB" + chTag + std::to_string(wh0Tag) + etaTag + IDTag])->Fill(probeMuTk.Eta(), dtPrim.bxTrackFinder());
				  
				  if (bx >= m_tnpConfig.probe_minPrimBX)
				    {
				      m_effs[TIMING]["earlyEffVsEtaMB" + chTag + std::to_string(wh0Tag) + etaTag + IDTag]->Fill(bx < 0, probeMuTk.Eta());
				      m_effs[TIMING]["lateEffVsEtaMB" + chTag + std::to_string(wh0Tag) + etaTag + IDTag]->Fill(bx > 0, probeMuTk.Eta());
				    }
				  
				  if ( ( wh0Tag==0 && 
					 probeMuTk.Eta() >  m_tnpConfig.probe_minEtaBX && 
					 probeMuTk.Eta() <  m_tnpConfig.probe_maxEtaBX
					 ) ||
				       ( wh0Tag==1 && 
					 probeMuTk.Eta() > -m_tnpConfig.probe_maxEtaBX && 
					 probeMuTk.Eta() < -m_tnpConfig.probe_minEtaBX
					 ) 
				       )
				    {

				      if(bx < firstBX[match.id_r - 1] &&
					 bx >= m_tnpConfig.probe_minPrimBX )
					firstBX[match.id_r - 1] = bx;
			      
				      if(probeMuTk.Pt() > 400.)
					{
					  static_cast<TH1F*>(m_histos[TIMING]["highPtTrigVsPhiMB" + chTag + etaTag + IDTag])->Fill(probeMuTk.Phi(), weight);
					  static_cast<TH1F*>(m_histos[TIMING]["highPtTrigBxMB" + chTag + etaTag + IDTag])->Fill(bx, weight);					  
					}
				      
				      static_cast<TProfile*>(m_histos[TIMING]["bxTrigVsPtMB" + chTag + etaTag + IDTag])->Fill(probeMuTk.Pt(), bx);
				      
				      if (bx >= m_tnpConfig.probe_minPrimBX)
					{ 
					  m_effs[TIMING]["bxm1EffVsPtMB" + chTag + etaTag + IDTag]->Fill(bx == -1, probeMuTk.Pt());
					  m_effs[TIMING]["bxm2EffVsPtMB" + chTag + etaTag + IDTag]->Fill(bx == -2, probeMuTk.Pt());
					  
					  m_effs[TIMING]["earlyEffVsPtMB" + chTag + etaTag + IDTag]->Fill(bx < 0, probeMuTk.Pt());
					  m_effs[TIMING]["earlyEffVsPhiMB" + chTag + etaTag + IDTag]->Fill(bx < 0, probeMuTk.Phi());

					  m_effs[TIMING]["lateEffVsPtMB" + chTag + etaTag + IDTag]->Fill(bx > 0, probeMuTk.Pt());
					  m_effs[TIMING]["lateEffVsPhiMB" + chTag + etaTag + IDTag]->Fill(bx > 0, probeMuTk.Phi());
					}
				      
				      m_histos[TIMING]["bxTrigMB" + chTag + etaTag + IDTag]->Fill(bx, weight);
				      m_histos[TRIG]["qualTrigMB" + chTag + etaTag + IDTag]->Fill(dtPrim.quality, weight);
				  
				    }
				}
			    }
			}

		      m_histos[TRIG]["nTrigMB1" + etaTag + IDTag]->Fill(nTrig[0], weight);
		      m_histos[TRIG]["nTrigMB2" + etaTag + IDTag]->Fill(nTrig[1], weight);
		      m_histos[TRIG]["nTrigMB3" + etaTag + IDTag]->Fill(nTrig[2], weight);
		      m_histos[TRIG]["nTrigMB4" + etaTag + IDTag]->Fill(nTrig[3], weight);
		      
		      if (firstBX[0] < 9)
			m_effs[TIMING]["firstEffVsPtMB1" + etaTag + IDTag]->Fill(firstBX[0] < 0,probeMuTk.Pt());
		      if (firstBX[1] < 9)
			m_effs[TIMING]["firstEffVsPtMB2" + etaTag + IDTag]->Fill(firstBX[1] < 0,probeMuTk.Pt());
		      if (firstBX[2] < 9)
			m_effs[TIMING]["firstEffVsPtMB3" + etaTag + IDTag]->Fill(firstBX[2] < 0,probeMuTk.Pt());
		      if (firstBX[3] < 9)
			m_effs[TIMING]["firstEffVsPtMB4" + etaTag + IDTag]->Fill(firstBX[3] < 0,probeMuTk.Pt());

		    }
		}
	    }
	}
    }
}


//Functions
void muon_pog::parseConfig(const std::string configFile, muon_pog::TagAndProbeConfig & tpConfig,
			   std::vector<muon_pog::SampleConfig> & sampleConfigs)
{

  boost::property_tree::ptree pt;
  
  try
    {
      boost::property_tree::ini_parser::read_ini(configFile, pt);
    }
  catch (boost::property_tree::ini_parser::ini_parser_error iniParseErr)
    {
      std::cout << "[TagAndProbeConfig] Can't open : " << iniParseErr.filename()
		<< "\n\tin line : " << iniParseErr.line()
		<< "\n\thas error :" << iniParseErr.message()
		<< std::endl;
      throw std::runtime_error("Bad INI parsing");
    }

  for( auto vt : pt )
    {
      if (vt.first.find("TagAndProbe") != std::string::npos)
	tpConfig = muon_pog::TagAndProbeConfig(vt);
      else
	sampleConfigs.push_back(muon_pog::SampleConfig(vt));
    }
}
