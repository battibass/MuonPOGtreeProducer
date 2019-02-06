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

    std::vector<TString> fileNames;  
    std::vector<double>  weights;
    std::vector<double>  pBins;
    std::vector<double>  weightBins;
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
    std::vector<int>     toArrayI(const std::string & entries); 
    std::vector<double>  toArrayF(const std::string & entries); 
    std::vector<TString> toArrayTS(const std::string & entries); 
    
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
    Float_t probe_minHighPt;      
    Int_t   probe_minNSeg;      

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

    enum HistoType { KIN=0, CONT, HIT, EFF, SEG};
      
    Plotter(muon_pog::TagAndProbeConfig tnpConfig, muon_pog::SampleConfig & sampleConfig) :
      m_tnpConfig(tnpConfig) , m_sampleConfig(sampleConfig) {};
    ~Plotter() {};
    
    void book(TFile *outFile);
    void fill(const std::vector<muon_pog::Muon> & muons, const muon_pog::HLT & hlt, const muon_pog::Event & ev, float weight);

    void fillBasicPlots(TString etaIdTag, TString categoryTag,
			const TLorentzVector & probeMuTk, const std::array<Int_t,8> & showers, double weight);

    void fillShowerPlots(TString etaIdTag, TString categoryTag, const muon_pog::Muon & muon,
			 const TLorentzVector & genMuTk, const std::array<bool,8> & showers, double weight);

    void bookBasic(TFile *outFile, const TString & etaIdTag,
		   const TString & sampleTag, const TString & categoryTag);

    void bookShower(TFile *outFile, const TString & etaIdTag,
		    const TString & sampleTag, const TString & categoryTag);

    std::map<Plotter::HistoType, std::map<TString, TH1* > > m_histos;
    std::map<Plotter::HistoType, std::map<TString, TEfficiency *> > m_effs;
    TagAndProbeConfig m_tnpConfig;
    SampleConfig m_sampleConfig;

  private:

    double weightFromP(double p);
        
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

      int iWeight = 0;
      for (auto fileName : plotter.m_sampleConfig.fileNames)
	{

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
	      float weight =plotter.m_sampleConfig.weights[iWeight];
	      
	      plotter.fill(ev->muons, ev->hlt, (*ev), weight);
	    }

	  delete ev;
	  inputFile->Close();
	  std::cout << std::endl;

	  ++iWeight;

	}
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

      probe_minNSeg   = vt.second.get<Int_t>("probe_minNSeg");
      probe_minHighPt = vt.second.get<Float_t>("probe_minHighPt");

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

std::vector<TString> muon_pog::TagAndProbeConfig::toArray(const std::string& entries)
{
  
  std::vector<TString> result;
  std::stringstream sentries(entries);
  std::string item;
  while(std::getline(sentries, item, ','))
    result.push_back(TString(item.c_str()));
  
  return result;
    
}
muon_pog::SampleConfig::SampleConfig(boost::property_tree::ptree::value_type & vt)
{
 
 try
    {
      fileNames = toArrayTS(vt.second.get<std::string>("fileNames"));
      weights = toArrayF(vt.second.get<std::string>("weights"));
      pBins = toArrayF(vt.second.get<std::string>("pBins"));
      weightBins = toArrayF(vt.second.get<std::string>("weightBins"));
      sampleName = TString(vt.first.c_str());
      cSection = vt.second.get<Float_t>("cSection");
      nEvents = vt.second.get<Float_t>("nEvents");
      applyReweighting = vt.second.get<Bool_t>("applyReweighting");
      runs = toArrayI(vt.second.get<std::string>("runs"));
    }
  
  catch (boost::property_tree::ptree_bad_data bd)
    {
      std::cout << "[TagAndProbeConfig] Can't get data : has error : "
		<< bd.what() << std::endl;
      throw std::runtime_error("Bad INI variables");
    }
  
}


std::vector<int> muon_pog::SampleConfig::toArrayI(const std::string& entries)
{
  
  std::vector<int> result;
  std::stringstream sentries(entries);
  std::string item;
  while(std::getline(sentries, item, ','))
    result.push_back(std::stoi(item.c_str()));
  
  return result;
    
}

std::vector<double> muon_pog::SampleConfig::toArrayF(const std::string& entries)
{
  
  std::vector<double> result;
  std::stringstream sentries(entries);
  std::string item;
  while(std::getline(sentries, item, ','))
    result.push_back(std::stod(item.c_str()));
  
  return result;
  
}

std::vector<TString> muon_pog::SampleConfig::toArrayTS(const std::string& entries)
{
  
  std::vector<TString> result;
  std::stringstream sentries(entries);
  std::string item;
  while(std::getline(sentries, item, ','))
    result.push_back(TString(item.c_str()));
  
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
  outFile->mkdir(sampleTag+"/hit_count");
  outFile->mkdir(sampleTag+"/control");
  outFile->mkdir(sampleTag+"/segment");
  
  for (auto fEtaBin : m_tnpConfig.probe_fEtaBins)
    {
      TString etaTag = "_fEtaMin" + fEtaBin.first + "_fEtaMax" + fEtaBin.second;
      
      for ( auto & probe_ID : m_tnpConfig.probe_IDs)
	{
	  TString IDTag = "_" + probe_ID;
	  
	  TString completeTag = etaTag + IDTag + sampleTag;	  
	  
	  outFile->cd(sampleTag+"/kinematical_variables");
	  
	  m_histos[KIN]["probeP" + etaTag + IDTag]  = new TH1F("probeP" + completeTag,
							       "probeP" + completeTag +
							       ";p (GeV);# entries",
							       250,0.,2500.);
	  
	  m_histos[KIN]["probeEta" + etaTag + IDTag] = new TH1F("probeEta" + completeTag,
								"probeEta" + completeTag +
								";#eta;# entries",
								48,-2.4, 2.4);
	  
	  m_histos[KIN]["probePhi" + etaTag + IDTag] = new TH1F("probePhi" + completeTag,
								"probePhi" + completeTag +
								";#phi;# entries",
								96,-TMath::Pi(),TMath::Pi());

 	  m_histos[KIN]["probePVsEta" + etaTag + IDTag]  = new TH2F("probePVsEta" + completeTag,
								    "probePVsEta" + completeTag +
								    "#eta ;p (GeV)",
								    12, 0., 2.4, 50,0.,2500.);

	  bookBasic(outFile, etaTag + IDTag, sampleTag, "Chamb");
	  bookBasic(outFile, etaTag + IDTag, sampleTag, "50");
	  bookBasic(outFile, etaTag + IDTag, sampleTag, "25");
	  bookBasic(outFile, etaTag + IDTag, sampleTag, "15");

	  bookBasic(outFile, etaTag + IDTag, sampleTag, "Segment");

	  bookShower(outFile, etaTag + IDTag, sampleTag, "Segment2");	  
	  bookShower(outFile, etaTag + IDTag, sampleTag, "Segment3");	  

	  std::vector<TString> analysisTags = { "1", "2", "3", "4", "5", "6", "7" };
	  
	  for ( const auto & analysisTag : analysisTags)
	    {
	      bookShower(outFile, etaTag + IDTag, sampleTag, "Chamb_" + analysisTag);	  
	      bookShower(outFile, etaTag + IDTag, sampleTag, "50_" + analysisTag);	  
	      bookShower(outFile, etaTag + IDTag, sampleTag, "25_" + analysisTag);	  
	      bookShower(outFile, etaTag + IDTag, sampleTag, "15_" + analysisTag);	  
	    }

	  std::vector<TString> chambTags = { "MB1", "MB2", "MB3", "MB4", "ME1", "ME2", "ME3", "ME4" };

	  outFile->cd(sampleTag+"/segment");

	  for ( const auto & chambTag : chambTags)
	    {
 	      m_histos[SEG]["errXVsPosX" + chambTag + etaTag + IDTag]  = new TH2F("errXVsPosX" + chambTag + etaTag + IDTag,
										 "errXVsPosX" + chambTag + etaTag + IDTag +
										 ";#Delta(pos_x); #Delta(err_x)",
										 201, -10.05, 10.05, 201, -1.005, 1.005);

 	      m_histos[SEG]["errXVsPhi" + chambTag + etaTag + IDTag]  = new TH2F("errXVsPhi" + chambTag + etaTag + IDTag,
										 "errXVsPhi" + chambTag + etaTag + IDTag +
										 ";#Delta(phi); #Delta(err_x)",
										 201, -.01005, .01005, 201, -1.005, 1.005);

	      m_histos[SEG]["posYVsPosX" + chambTag + etaTag + IDTag]  = new TH2F("posYVsPosX" + chambTag + etaTag + IDTag,
										  "posYVsPosX" + chambTag + etaTag + IDTag +
										  ";#Delta(pos_x); #Delta(pos_y)",
										  201, -100.5, 100.5, 201, -100.5, 100.5);

	      m_histos[SEG]["nHitsVsPosX" + chambTag + etaTag + IDTag]  = new TH2F("nHitsVsPosX" + chambTag + etaTag + IDTag,
										  "nHitsVsPosX" + chambTag + etaTag + IDTag +
										  ";#Delta(pos_x); #Delta(# hits)",
										  201, -10.05, 10.05, 21, -10.5, 10.5);

	      m_histos[SEG]["nDigiVsNSegShower" + chambTag + etaTag + IDTag]  = new TH2F("nDigiVsNSegShower" + chambTag + etaTag + IDTag,
											 "nDigiVsNSegShower" + chambTag + etaTag + IDTag +
											 ";# segment ;# digi",
											 101, -0.5, 100.5, 101, -0.5, 100.5);

	      m_histos[SEG]["nDigiVsNSegAll" + chambTag + etaTag + IDTag]  = new TH2F("nDigiVsNSegAll" + chambTag + etaTag + IDTag,
										      "nDigiVsNSegAll" + chambTag + etaTag + IDTag +
										      ";# segment ;# digi",
										      101, -0.5, 100.5, 101, -0.5, 100.5);
	    }

	}

    }

  for ( auto & pair : m_effs[EFF])
    {
      pair.second->SetUseWeightedEvents();
      pair.second->SetStatisticOption(TEfficiency::kFNormal);
    }  

}

void muon_pog::Plotter::bookBasic(TFile *outFile,
				  const TString & etaIdTag,
				  const TString & sampleTag,
				  const TString & categoryTag)
{

  outFile->cd(sampleTag+"/hit_count");

  // const Double_t ptBins[12] = {0., 50., 100., 150., 200., 300., 400., 500., 600., 700., 800., 1200.};
  std::vector<std::string> chambers = { "MB1", "MB2", "MB3", "MB4", "ME1", "ME2", "ME3", "ME4" };
  for (const auto & chamb : chambers)
    {

      TString plotTag = chamb + categoryTag + etaIdTag;
      TString completeTag = chamb + categoryTag + etaIdTag + sampleTag;

      m_histos[HIT]["nHit" + plotTag] = new TH1F("nHit" + completeTag, 
						 "nHit" + completeTag +
						 ";# hits per station; # entries", 
						  200, -0.5, 199.5);

      m_histos[HIT]["nHitVsP" + plotTag] = new TProfile("nHitVsP" + completeTag, 
							 "nHitVsP" + completeTag +
							 ";muon p (GeV/c);# hits per station", 
							 12, 0., 2400., -0.5, 199.5);

      m_histos[HIT]["nHitVsPhi" + plotTag] = new TProfile("nHitVsPhi" + completeTag, 
							  "nHitVsPhi" + completeTag +
							  ";muon #phi;# hits per station", 
							  96, -TMath::Pi(), TMath::Pi(), -0.5, 199.5);

      m_histos[HIT]["nHitVsEta" + plotTag] = new TProfile("nHitVsEta" + completeTag, 
							  "nHitVsEta" + completeTag +
							  ";muon #eta;# hits per station", 
							  48, -2.4, 2.4, -0.5, 199.5);
 
    }

  TString plotTag = categoryTag + etaIdTag;
  TString completeTag = categoryTag + etaIdTag + sampleTag;

  m_histos[HIT]["nDtHitVsP" + plotTag] = new TProfile("nDtHitsVsP" + completeTag, 
						      "nDtHitsVsP" + completeTag +
						      ";muon p (GeV/c);hits per station", 
						      24, 0., 2400., -0.5, 199.5);
  
  m_histos[HIT]["nCscHitVsP" + plotTag] = new TProfile("nCscHitsVsP" + completeTag, 
						       "nCscHitsVsP" + completeTag +
						       ";muon p (GeV/c);hits per station", 
						       24, 0., 2400., -0.5, 199.5);
  
}


void muon_pog::Plotter::bookShower(TFile *outFile,
				   const TString & etaIdTag,
				   const TString & sampleTag,
				   const TString & categoryTag)
{

  outFile->cd(sampleTag+"/efficiencies");

  // const Double_t ptBins[12] = {0., 50., 100., 150., 200., 300., 400., 500., 600., 700., 800., 1200.};
  std::vector<std::string> chambers = { "MB1", "MB2", "MB3", "MB4", "ME1", "ME2", "ME3", "ME4" };
  for (const auto & chamb : chambers)
    {

      TString plotTag = chamb + categoryTag + etaIdTag;
      TString completeTag = chamb + categoryTag + etaIdTag + sampleTag;

      m_effs[EFF]["nShowersVsP" + plotTag] = new TEfficiency("nShowersVsP" + completeTag, 
							      "nShowersVsP" + completeTag +
							      ";muon p (GeV/c);shower probability", 
							      12, 0., 2400.);

      m_effs[EFF]["nShowersVsPhi" + plotTag] = new TEfficiency("nShowersVsPhi" + completeTag, 
							       "nShowersVsPhi" + completeTag +
							       ";muon #phi;shower probability", 
							       96, -TMath::Pi(), TMath::Pi());

      m_effs[EFF]["nShowersVsEta" + plotTag] = new TEfficiency("nShowersVsEta" + completeTag, 
							       "nShowersVsEta" + completeTag +
							       ";muon #eta;shower probability", 
							       48, -2.4, 2.4);
 
    }

  TString plotTag = categoryTag + etaIdTag;
  TString completeTag = categoryTag + etaIdTag + sampleTag;

  m_effs[EFF]["nDtShowersVsP" + plotTag] = new TEfficiency("nDtShowersVsP" + completeTag, 
							   "nDtShowersVsP" + completeTag +
							   ";muon p (GeV/c);shower probability", 
							   24, 0., 2400.);
  
  m_effs[EFF]["nDtShowersVsPhi" + plotTag] = new TEfficiency("nDtShowersVsPhi" + completeTag, 
							     "nDtShowersVsPhi" + completeTag +
							     ";muon #phi;shower probability", 
							     96, -TMath::Pi(), TMath::Pi());
  
  m_effs[EFF]["nDtShowersVsEta" + plotTag] = new TEfficiency("nDtShowersVsEta" + completeTag, 
							     "nDtShowersVsEta" + completeTag +
							     ";muon #eta;shower probability", 
							     48, -2.4, 2.4);
  
  m_effs[EFF]["nDtShowersPerCh" + plotTag] = new TEfficiency("nDtShowersPerCh" + completeTag, 
							     "nDtShowersPerCh" + completeTag +
							     ";station # ;shower probability", 
							     4, 0.5, 4.5);

  m_effs[EFF]["nCscShowersVsP" + plotTag] = new TEfficiency("nCscShowersVsP" + completeTag, 
							    "nCscShowersVsP" + completeTag +
							    ";muon p (GeV/c);shower probability", 
							    24, 0., 2400.);
  
  m_effs[EFF]["nCscShowersVsPhi" + plotTag] = new TEfficiency("nCscShowersVsPhi" + completeTag, 
							      "nCscShowersVsPhi" + completeTag +
							      ";muon #phi;shower probability", 
							      96, -TMath::Pi(), TMath::Pi());
  
  m_effs[EFF]["nCscShowersVsEta" + plotTag] = new TEfficiency("nCscShowersVsEta" + completeTag, 
							      "nCscShowersVsEta" + completeTag +
							      ";muon #eta;shower probability", 
							      48, -2.4, 2.4);
  
  m_effs[EFF]["nCscShowersPerCh" + plotTag] = new TEfficiency("nCscShowersPerCh" + completeTag, 
							      "nCscShowersPerCh" + completeTag +
							      ";station # ;shower probability", 
							      4, 0.5, 4.5);

  outFile->cd(sampleTag+"/control");
  
  m_histos[CONT]["nDtShowers" + plotTag] = new TH1F("nDtShowers" + completeTag,
						    "nDtShowers" + completeTag +
						    ";# of stations with showers;# entries",
						    4, 0.5, 4.5);
  
  m_histos[CONT]["nCscShowers" + plotTag] = new TH1F("nCscShowers" + completeTag,
						     "nCscShowers" + completeTag +
						     ";# of stations with showers;# entries",
						     4, 0.5, 4.5);

  m_histos[CONT]["nDtMatchStVsShowers" + plotTag] = new TProfile("nDtMatchStVsShowers" + completeTag,
								 "nDtMatchStVsShowers" + completeTag +
								 ";# of stations with showers;# of matched stations",
								 4, -0.5, 3.5, -0.5, 4.5);
  
  m_histos[CONT]["nCscMatchStVsShowers" + plotTag] = new TProfile("nCscMatchStVsShowers" + completeTag,
								  "nCscMatchStVsShowers" + completeTag +
								  ";# of stations with showers;# of matched stations",
								  4, -0.5, 3.5, -0.5, 4.5);
  
}


void muon_pog::Plotter::fill(const std::vector<muon_pog::Muon> & muons,
			     const muon_pog::HLT & hlt, const muon_pog::Event & ev, float weight)
{

  bool isGoodRun = false;

  if (m_sampleConfig.runs.at(0) == 0 ||
      (ev.runNumber >= m_sampleConfig.runs.at(0) &&
       ev.runNumber <  m_sampleConfig.runs.at(1) )
      )
    {
      isGoodRun = true;
    }

  if (!isGoodRun) return;
  
  std::vector<const muon_pog::Muon *> probeMuons;
  
  for (auto & muon : muons)
    {
      
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
	 muon.trackerIso < m_tnpConfig.probe_isoCutAbs &&
	 muon.trackerIso / muTk.Pt() < m_tnpConfig.probe_isoCut)
	{
	  probeMuons.push_back(&muon);
	}
	 
    }
  
  for (auto probeMuonPointer : probeMuons)
    {
      const muon_pog::Muon & probeMuon = *probeMuonPointer;
      
      TLorentzVector probeMuTk(muonTk(probeMuon,m_tnpConfig.muon_trackType));

      float weight = weightFromP(probeMuTk.P());

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
		      m_histos[KIN]["probeP" + etaTag + IDTag]->Fill(probeMuTk.P(), weight);
		      m_histos[KIN]["probeEta" + etaTag + IDTag]->Fill(probeMuTk.Eta(), weight);
		      m_histos[KIN]["probePhi" + etaTag + IDTag]->Fill(probeMuTk.Phi(), weight);
		      static_cast<TH2F*>(m_histos[KIN]["probePVsEta" + etaTag + IDTag])->Fill(probeMuTk.Eta(), std::abs(probeMuTk.P()), weight);

		      auto nShowersSegment = showerPerCh(probeMuon, ev.dtSegments, ev.cscSegments, 25.);

		      auto nShowersChamb = showerPerCh(probeMuon, 0);
		      auto nShowers50    = showerPerCh(probeMuon, 1);
		      auto nShowers25    = showerPerCh(probeMuon, 2);
		      auto nShowers15    = showerPerCh(probeMuon, 3);

		      fillBasicPlots(etaTag + IDTag, TString("Segment"), probeMuTk, nShowersSegment, weight);

		      fillBasicPlots(etaTag + IDTag, TString("Chamb"), probeMuTk, nShowersChamb, weight);
		      fillBasicPlots(etaTag + IDTag, TString("50"), probeMuTk, nShowers50, weight);
		      fillBasicPlots(etaTag + IDTag, TString("25"), probeMuTk, nShowers25, weight);
		      fillBasicPlots(etaTag + IDTag, TString("15"), probeMuTk, nShowers15, weight);

		      std::vector<TString> analysisTags = { "1", "2", "3", "4", "5", "6", "7" };
		      std::vector<int> cscCutTags = { 27, 36, 45, 54, 63, 72, 81 };
		      std::vector<int> dtCutTags  = { 12, 16, 20, 24, 28, 32, 36 };

		      auto analysisTagIt  = analysisTags.begin();
		      auto analysisTagEnd = analysisTags.end();

		      auto cscCutTagIt  = cscCutTags.begin();
		      auto cscCutTagEnd = cscCutTags.end();

		      auto dtCutTagIt  = dtCutTags.begin();
		      auto dtCutTagEnd = dtCutTags.end();

		      auto hasShowersSegment2 = hasShowerPerCh(probeMuon, ev.dtSegments, ev.cscSegments, 25., 2, 2);
		      auto hasShowersSegment3 = hasShowerPerCh(probeMuon, ev.dtSegments, ev.cscSegments, 25., 3, 3);

		      fillShowerPlots(etaTag + IDTag, TString("Segment2"), probeMuon, probeMuTk, hasShowersSegment2, weight);
		      fillShowerPlots(etaTag + IDTag, TString("Segment3"), probeMuon, probeMuTk, hasShowersSegment3, weight);

		      for (; analysisTagIt != analysisTagEnd && 
			     cscCutTagIt   != cscCutTagEnd && 
			     dtCutTagIt    != dtCutTagEnd; 
			     ++analysisTagIt, ++dtCutTagIt, ++cscCutTagIt)
			{

			  auto hasShowersChamb = hasShowerPerCh(probeMuon, (*dtCutTagIt), (*cscCutTagIt), 0);
			  auto hasShowers50    = hasShowerPerCh(probeMuon, (*dtCutTagIt), (*cscCutTagIt), 1);
			  auto hasShowers25    = hasShowerPerCh(probeMuon, (*dtCutTagIt), (*cscCutTagIt), 2);
			  auto hasShowers15    = hasShowerPerCh(probeMuon, (*dtCutTagIt), (*cscCutTagIt), 3);

			  fillShowerPlots(etaTag + IDTag, TString("Chamb_") + (*analysisTagIt), probeMuon, probeMuTk, hasShowersChamb, weight);
			  fillShowerPlots(etaTag + IDTag, TString("50_") + (*analysisTagIt), probeMuon, probeMuTk, hasShowers50, weight);		      
			  fillShowerPlots(etaTag + IDTag, TString("25_") + (*analysisTagIt), probeMuon, probeMuTk, hasShowers25, weight);		      
			  fillShowerPlots(etaTag + IDTag, TString("15_") + (*analysisTagIt), probeMuon, probeMuTk, hasShowers15, weight);

			}

		      auto segmentsPerCh  = getSegmentsPerCh(probeMuon, ev.dtSegments, ev.cscSegments);
		
		      Int_t iCh = 0;
		      std::vector<TString> chambTags = { "MB1", "MB2", "MB3", "MB4", "ME1", "ME2", "ME3", "ME4" };

		      for ( const auto & chambTag : chambTags)
			{

			  const auto & segments = segmentsPerCh[iCh]; 
			  
			  if(nShowers25[iCh] > 0)
			    {
			      m_histos[SEG]["nDigiVsNSegShower" + chambTag + etaTag + IDTag]->Fill(nShowers25[iCh],nShowersSegment[iCh]);
			      m_histos[SEG]["nDigiVsNSegAll" + chambTag + etaTag + IDTag]->Fill(nShowers25[iCh],segments.size());
			    }

			  auto seg1It = segments.cbegin();
			  auto segEnd = segments.cend();
			  
			  // std::cout << chambTag << std::endl;

			  for (; seg1It != segEnd; ++seg1It)
			    {
			      auto seg2It = seg1It;
			      seg2It++;
			      for (; seg2It != segEnd; ++seg2It)
				{				      
				  Float_t dX    = seg1It->x - seg2It->x;
				  Float_t dY    = seg1It->y - seg2It->y;
				  Float_t dPhi  = seg1It->phi - seg2It->phi;
				  Float_t dErrX = seg1It->errx - seg2It->errx;
				  Float_t dNHit = seg1It->nHitsX - seg2It->nHitsX;

				  // std::cout << (seg1It == seg2It) << " " << dX << " " << dY << " " << dErrX << " " << dNHit << std::endl;
	  
				  m_histos[SEG]["errXVsPosX" + chambTag + etaTag + IDTag]->Fill(dX,dErrX);

				  m_histos[SEG]["errXVsPhi" + chambTag + etaTag + IDTag]->Fill(dPhi,dErrX);
				  
				  m_histos[SEG]["posYVsPosX" + chambTag + etaTag + IDTag]->Fill(dX,dY);

				  m_histos[SEG]["nHitsVsPosX" + chambTag + etaTag + IDTag]->Fill(dX,dNHit);
				}
			    }

			  ++iCh;

			  // std::cout << std::endl << std::endl;

			}
		      
		    }
		}
	    }
	}
    }
}

void muon_pog::Plotter::fillBasicPlots(TString etaIdTag,
				       TString categoryTag,
				       const TLorentzVector & probeMuTk,
				       const std::array<Int_t,8> & showers,
				       double weight)
{

  Int_t iChamb = 0;
  std::vector<std::string> chambers = { "MB1", "MB2", "MB3", "MB4", "ME1", "ME2", "ME3", "ME4" };
  for (const auto & chamb : chambers)
    {

      if (showers[iChamb] > 0)
	{
	  m_histos[HIT]["nHit" + chamb + categoryTag + etaIdTag]->Fill(showers[iChamb], weight);
	  static_cast<TProfile*>(m_histos[HIT]["nHitVsP" + chamb + categoryTag + etaIdTag])->Fill(probeMuTk.P(), showers[iChamb], weight);
	  static_cast<TProfile*>(m_histos[HIT]["nHitVsPhi" + chamb + categoryTag + etaIdTag])->Fill(probeMuTk.Phi(), showers[iChamb], weight);
	  static_cast<TProfile*>(m_histos[HIT]["nHitVsEta" + chamb + categoryTag + etaIdTag])->Fill(probeMuTk.Eta(), showers[iChamb], weight);
	  if (std::abs(probeMuTk.Eta()) < 0.9)
	    static_cast<TProfile*>(m_histos[HIT]["nDtHitVsP" + categoryTag + etaIdTag])->Fill(probeMuTk.P(),showers[iChamb], weight);
	  if (std::abs(probeMuTk.Eta()) > 1.2 &&
	      std::abs(probeMuTk.Eta()) < 2.4)
	    static_cast<TProfile*>(m_histos[HIT]["nCscHitVsP" + categoryTag + etaIdTag])->Fill(probeMuTk.P(),showers[iChamb], weight);
	}

      iChamb++;
    }

}

void muon_pog::Plotter::fillShowerPlots(TString etaIdTag, 
					TString categoryTag,
					const muon_pog::Muon & muon,
					const TLorentzVector & refMuTk,
					const std::array<bool,8> & showers,
					double weight)
{


  if (std::abs(refMuTk.Eta()) < 0.9) {

    Int_t iChamb = 0;
    std::vector<std::string> chambers = { "MB1", "MB2", "MB3", "MB4"};
    for (const auto & chamb : chambers)
      {
	m_effs[EFF]["nShowersVsP" + chamb + categoryTag + etaIdTag]->FillWeighted(showers[iChamb], weight, refMuTk.P());
	m_effs[EFF]["nShowersVsPhi" + chamb + categoryTag + etaIdTag]->FillWeighted(showers[iChamb], weight, refMuTk.Phi());
	m_effs[EFF]["nShowersVsEta" + chamb + categoryTag + etaIdTag]->FillWeighted(showers[iChamb], weight, refMuTk.Eta());
	
	m_effs[EFF]["nDtShowersPerCh" + categoryTag + etaIdTag]->FillWeighted(showers[iChamb], weight, iChamb + 1);

	iChamb++;
      }

    int nDtShowers = 0 +
      (showers[0] ? 1 : 0) +
      (showers[1] ? 1 : 0) +
      (showers[2] ? 1 : 0) +
      (showers[3] ? 1 : 0);
  
    m_effs[EFF]["nDtShowersVsP" + categoryTag + etaIdTag]->FillWeighted(nDtShowers, weight, refMuTk.P());
    m_effs[EFF]["nDtShowersVsPhi" + categoryTag + etaIdTag]->FillWeighted(nDtShowers, weight, refMuTk.Phi());
    m_effs[EFF]["nDtShowersVsEta" + categoryTag + etaIdTag]->FillWeighted(nDtShowers, weight, refMuTk.Eta());
    
    m_histos[CONT]["nDtShowers" + categoryTag + etaIdTag]->Fill(nDtShowers, weight);
    
    static_cast<TProfile*>(m_histos[CONT]["nDtMatchStVsShowers" + categoryTag + etaIdTag])->Fill(nDtShowers,
												 muon.trkMuonMatchedStations, 
												 weight);
  }


  if (std::abs(refMuTk.Eta()) > 1.2 &&
      std::abs(refMuTk.Eta()) < 2.4) {

    Int_t iChamb = 4;
    std::vector<std::string> chambers = { "ME1", "ME2", "ME3", "ME4" };
    for (const auto & chamb : chambers)
      {
	m_effs[EFF]["nShowersVsP" + chamb + categoryTag + etaIdTag]->FillWeighted(showers[iChamb], weight, refMuTk.P());
	m_effs[EFF]["nShowersVsPhi" + chamb + categoryTag + etaIdTag]->FillWeighted(showers[iChamb], weight, refMuTk.Phi());
	m_effs[EFF]["nShowersVsEta" + chamb + categoryTag + etaIdTag]->FillWeighted(showers[iChamb], weight, refMuTk.Eta());
	
	m_effs[EFF]["nCscShowersPerCh" + categoryTag + etaIdTag]->FillWeighted(showers[iChamb], weight, iChamb - 3);

	iChamb++;
      }

    int nCscShowers = 0 +
      (showers[4] ? 1 : 0) +
      (showers[5] ? 1 : 0) +
      (showers[6] ? 1 : 0) +
      (showers[7] ? 1 : 0);
    
    m_effs[EFF]["nCscShowersVsP" + categoryTag + etaIdTag]->FillWeighted(nCscShowers, weight, refMuTk.P());
    m_effs[EFF]["nCscShowersVsPhi" + categoryTag + etaIdTag]->FillWeighted(nCscShowers, weight, refMuTk.Phi());
    m_effs[EFF]["nCscShowersVsEta" + categoryTag + etaIdTag]->FillWeighted(nCscShowers, weight, refMuTk.Eta());

    m_histos[CONT]["nCscShowers" + categoryTag + etaIdTag]->Fill(nCscShowers, weight);

    static_cast<TProfile*>(m_histos[CONT]["nCscMatchStVsShowers" + categoryTag + etaIdTag])->Fill(nCscShowers,
												  muon.trkMuonMatchedStations, 
												  weight);
  }

}


double muon_pog::Plotter::weightFromP(double p)
{
  double weight = 0;

  auto pBinIt  = m_sampleConfig.pBins.begin();
  auto pBinEnd = m_sampleConfig.pBins.end();

  auto weightIt  = m_sampleConfig.weightBins.begin();
  auto weightEnd = m_sampleConfig.weightBins.end();

  for (; pBinIt != pBinEnd && weightIt != weightEnd; ++pBinIt, ++ weightIt)
    {
      if (p <= (*pBinIt))
	break;
      weight = (*weightIt);
    }

  return weight;

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
