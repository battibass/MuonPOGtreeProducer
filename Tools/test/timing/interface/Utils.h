#ifndef Utils_h__
#define Utils_h__

#include "MuonPogTree.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TChain.h"

#include <iostream>
#include <string>

#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>

namespace muon_pog
{

  // The names says it all
  double deltaR(double eta1, double phi1, double eta2, double phi2)
  {
    double deta = eta1 - eta2;
    double dphi = phi1 - phi2;
    while (dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
    while (dphi <= -TMath::Pi()) dphi += 2*TMath::Pi();

    return sqrt(deta*deta + dphi*dphi);
  }

    
  // Check if a muon_pog::Muon passes a given ID (embedded using official selectors)
  // Valid IDs are GLOBAL, TRACKER, SOFT, LOOSE, MEDIUM, TIGHT, HIGHPT
  bool hasGoodId(const muon_pog::Muon & muon, TString muId)
  {
    
    if (muId == "GLOBAL")      return muon.isGlobal  == 1;
    else if (muId == "TRACKER")return muon.isTracker == 1;
    else if (muId == "TIGHT")  return muon.isTight   == 1;
    else if (muId == "MEDIUM") return muon.isMedium  == 1;
    else if (muId == "LOOSE")  return muon.isLoose   == 1;
    else if (muId == "HIGHPTPOG") return muon.isHighPt  == 1;
    else if (muId == "HIGHPT") return muon.isGlobal &&
				            muon.isTracker &&
				            muon.glbMuonValidHits  > 0 &&
				            muon.trkPixelValidHits > 0 &&
				            muon.trkTrackerLayersWithMeas > 5 &&
				            muon.trkMuonZPrimeMatchedStations == 1 &&
				            std::abs(muon.dxy) < 0.2 &&
				            std::abs(muon.dz)  < 0.5 &&
				            muon.fitPtErr(muon_pog::MuonFitType::TUNEP) /
				            muon.fitPt(muon_pog::MuonFitType::TUNEP) < 0.3;
    else if (muId == "SOFT")   return muon.isSoft == 1;
    else
      {
	std::cout << "[Plotter::hasGoodId]: Invalid muon id : "
		  << muId << std::endl;
	exit(900);
      }
    
    return 0;
 
  };

  std::array<Int_t,4> showersPerCh(const muon_pog::Muon & muon,
				   const std::vector<muon_pog::MuonSegment> & dtSegments,
				   Float_t deltaPhi)
    {
        
      std::array<Int_t,4> nExtraSegPerCh = {0, 0, 0, 0};
      std::map<Int_t,std::vector<muon_pog::MuonSegment> > dtSegmentPerCh;
        
        
      for (const auto & match : muon.matches)
        {
	  
	  Int_t ch = match.id_r;
            
	  if (match.type != muon_pog::MuonDetType::DT)
	    continue;
            
	  std::vector<std::size_t>::const_iterator indexIt  = match.indexes.begin();
	  std::vector<std::size_t>::const_iterator indexEnd = match.indexes.end();
        
	  std::vector<std::bitset<4> >::const_iterator qualIt  = match.matchQuals.begin();
	  std::vector<std::bitset<4> >::const_iterator qualEnd = match.matchQuals.end();
            
	  for (; qualIt != qualEnd && indexIt != indexEnd; ++indexIt, ++ qualIt)
            {
	      
	      std::bitset<4> mask(std::string("0010"));
	      std::bitset<4> mask1(std::string("0010"));

		if( !(mask & (*qualIt)).count() &&
		  std::abs(match.phi - dtSegments.at(*indexIt).phi) < deltaPhi)
		dtSegmentPerCh[ch - 1].push_back(dtSegments.at((*indexIt)));
            }

        }
     
      for (auto & pair : dtSegmentPerCh)
	{ 

      	  std::vector<muon_pog::MuonSegment>::iterator seg1It  = pair.second.begin();
            
      	  for (; seg1It != pair.second.end(); ++seg1It)
       	    {
       	      std::vector<muon_pog::MuonSegment>::iterator seg2It  = seg1It;
       	      seg2It++;
                    
       	      for (; seg2It != pair.second.end(); ++seg2It)
       		{
		  
       		  if(seg1It->id_eta == seg2It->id_eta &&
       		     seg1It->id_phi == seg2It->id_phi &&
       		     seg1It->id_r   == seg2It->id_r   &&
       		     ( std::abs(seg1It->x - seg2It->x)       < 0.01 ||
       		       std::abs(seg1It->errx - seg2It->errx) < 0.01
       		       ) &&
       		     seg1It->nHitsX == seg2It->nHitsX)
       		    {
       		      auto copy = seg2It;
		      pair.second.erase(copy);
		      --seg2It;
       		    }
       		}
       	    }
 
	}
      
      for (const auto & pair : dtSegmentPerCh)
        {
	  nExtraSegPerCh[pair.first] = pair.second.size();
        }

      return nExtraSegPerCh;

    };

  std::array<bool,4> hasShowerPerCh(const muon_pog::Muon & muon,
				    const std::vector<muon_pog::MuonSegment> & dtSegments,
				    const std::vector<muon_pog::DtDigiSummary> & dtDigis;
				    Float_t deltaPhi, Int_t nSeg)
    {
        
      std::array<bool,4> showerPerCh = {false, false, false, false};
      
      auto nExtraSegPerCh = showersPerCh(muon,dtSegments,deltaPhi)
	
	for (const auto & match : muon.matches)
	  {
	    
	    Int_t ch = match.id_r;
            
	    if (match.type != muon_pog::MuonDetType::DT)
	      continue;
	    
	    if (nExtraSegPerCh[ch - 1] >= nSeg)
	      {
		showersPerCh[ch - 1] = true;
		continue;
	      }
	    
	    for (const auto & digiSummary : dtDigis)
	      {
		if(digiSummary.id_phi == match.id_phi &&
		   digiSummary.id_eta == match.id_eta &&
		   digiSummary.id_r   == match.id_r   &&
		   ( digiSummary.n_phi1 >= 50 ||
		     digiSummary.n_phi2 >= 50 )
		   )
		  {
		    showersPerCh[ch - 1] = true;
		    continue;
		  }
	      }
	    
	  }
      
      return showerPerCh;
      
    };

  // Returns the charge muon_pog::Muon for a given fit 
  // Valid track fits are: PF, TUNEP, GLB, INNER, PICKY, DYT, TPFMS
  Int_t chargeFromTrk(const muon_pog::Muon & muon, 
		      const std::string & trackType)
  {

    if (trackType == "PF")         return muon.fits.at(muon_pog::MuonFitType::DEFAULT).charge;
    else if (trackType == "TUNEP") return muon.fits.at(muon_pog::MuonFitType::TUNEP).charge;
    else if (trackType == "GLB")   return muon.fits.at(muon_pog::MuonFitType::GLB).charge;
    else if (trackType == "INNER") return muon.fits.at(muon_pog::MuonFitType::INNER).charge;
    else if (trackType == "PICKY") return muon.fits.at(muon_pog::MuonFitType::PICKY).charge;
    else if (trackType == "DYT") return muon.fits.at(muon_pog::MuonFitType::DYT).charge;
    else if (trackType == "TPFMS") return muon.fits.at(muon_pog::MuonFitType::TPFMS).charge;
      
    else
      {
	std::cout << "[Plotter::chargeFromTrk]: Invalid track type: "
		  << trackType << std::endl;
	exit(900);
      }

    return 999;
    
  }


  // Return a TLorentz vector out of a given fit from muon_pog::Muon 
  // Valid track fits are: PF, TUNEP, GLB, INNER, PICKY, DYT, TPFMS
  TLorentzVector muonTk(const muon_pog::Muon & muon, 
			const std::string & trackType)
  {

    TLorentzVector result;
    muon_pog::MuonFit fit;
    if (trackType == "PF")
      fit = muon.fits.at(muon_pog::MuonFitType::DEFAULT);
    else if (trackType == "TUNEP")
      fit = muon.fits.at(muon_pog::MuonFitType::TUNEP);
    else if (trackType == "GLB")
      fit = muon.fits.at(muon_pog::MuonFitType::GLB);
    else if (trackType == "INNER")
      fit = muon.fits.at(muon_pog::MuonFitType::INNER);
    else if (trackType == "PICKY")
      fit = muon.fits.at(muon_pog::MuonFitType::PICKY);
    else if (trackType == "DYT")
      fit = muon.fits.at(muon_pog::MuonFitType::DYT);
    else if (trackType == "TPFMS")
      fit = muon.fits.at(muon_pog::MuonFitType::TPFMS);
    else
      {
	std::cout << "[Plotter::muonTk]: Invalid track type: "
		  << trackType << std::endl;
	exit(900);
      }

    result.SetPtEtaPhiM(fit.pt,fit.eta,fit.phi,.10565);
    return result;
    
  }


  // Checks if a trigger path fired using muon_pog::HLT
  // if the path name is "none" returns always true
  bool pathHasFired(const muon_pog::HLT  & hlt, std::string pathName)
  {
    
    if (pathName == "none")
      return true;
    
    for (auto path : hlt.triggers)
      {
	if (path.find(pathName) != std::string::npos)
	  {
	    return true;
	  }
      }
  }
  

  // Checks if the iner track of muon_pog::Muon matches geometrically 
  // in dR with a given muon_pog::HLT object filter
  // if the filter name is "none" returns always true
  bool hasFilterMatch(const muon_pog::Muon & muon,
		      const muon_pog::HLT  & hlt,
		      std::string & filter, Float_t dR)
  {

    if (filter == "none")
      return true;
    
    if (!hasGoodId(muon,"GLOBAL") &&
	!hasGoodId(muon,"TRACKER") ) 
      return false;

    TLorentzVector muTk = muonTk(muon,std::string("INNER"));

    for (auto object : hlt.objects)
      {
	if (object.filterTag.find(filter) != std::string::npos &&
	    deltaR(muTk.Eta(), muTk.Phi(), object.eta, object.phi) < dR)
	  return true;
      }
    
    return false;
    
  }

  bool hasMother(const muon_pog::GenParticle & gen, Int_t pdgId)
    {
      
      for (auto motherId : gen.mothers)
	{
	  if (abs(motherId) == pdgId)
	    {
	      return true;
	    }
	}
      
      return false;
      
    }


  const muon_pog::GenParticle * hasGenMatch(const muon_pog::Muon & muon,
					    const std::vector<muon_pog::GenParticle>  & gens,
					    Float_t dRCut, Int_t motherPdgId = 0, Int_t vetoPdgId = 0)
  {
   
    TLorentzVector muTk = muonTk(muon,std::string("INNER"));
    
    const muon_pog::GenParticle * bestGen = 0;
    Float_t bestDr = 999.;
    
    for (auto & gen : gens)
      {
	if (fabs(gen.pdgId) == 13)
	  {
	    Float_t dr = deltaR(muTk.Eta(), muTk.Phi(), gen.eta, gen.phi);
	    if (dr < dRCut && dr < bestDr && 
		((vetoPdgId == 0) || !hasMother(gen,vetoPdgId)))
	      {
		bestGen = &gen;
		bestDr = dr;
	      }
	  }
      }

    return (bestGen && ((motherPdgId == 0) || hasMother(*bestGen,motherPdgId))) ? bestGen : 0;
  }


  std::string exec(const char* cmd) {
    char buffer[100000];
    std::string result = "";
    std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");
    while (!feof(pipe.get())) {
        if (fgets(buffer, 100000, pipe.get()) != NULL)
            result += buffer;
    }
    return result;
  }
  
  TChain * openFileOrDir(std::string path) {
     
     bool isDir = false;
     if(path.find(".root") == std::string::npos) isDir = true;
     bool isEOS = false;
     if (path.find("/eos/cms/store/") != std::string::npos) isEOS = true;

     TChain *muonChain = new TChain("MuonPogTree/MUONPOGTREE", "MUONPOGTREE");
     std::vector<std::string> files;
     if(isDir && !isEOS) {
       std::string output = exec(("ls " + path + " | grep .root").c_str());
       std::string suboutput = output;
       while(suboutput.find(".root") != std::string::npos) {
          std::string theoutput = suboutput.substr(0, suboutput.find(".root")+5);
          suboutput = suboutput.substr(suboutput.find(".root") + 6, suboutput.size());
          if(path.at(path.size()-1) != '/') path = path + "/"; 
          std::string theeosoutput = path + theoutput;
          std::cout << "Adding file " << theeosoutput << std::endl;
          files.push_back(theeosoutput); 
       } 
     } else if (isDir && isEOS) {
       std::string output = exec(("/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select ls " + path + " | grep .root").c_str());
       std::string suboutput = output;
       while(suboutput.find(".root") != std::string::npos) {
          std::string theoutput = suboutput.substr(0, suboutput.find(".root")+5);
          suboutput = suboutput.substr(suboutput.find(".root") + 6, suboutput.size());
          if(path.at(path.size()-1) != '/') path = path + "/"; 
          std::string theeosoutput = "root://eoscms/" + path + theoutput;
          std::cout << "Adding file " << theeosoutput << std::endl;
          files.push_back(theeosoutput); 
       }
     } else if (!isDir && !isEOS) {
          files.push_back(path);
          std::cout << "Adding file " << path << std::endl;
     } else {
          files.push_back("root://eoscms/" + path);
          std::cout << "Adding file " << "root://eoscms/" << path << std::endl;
     }
     for(size_t i = 0; i < files.size(); i++) muonChain->Add(files[i].c_str());

     return muonChain;
  }
  
}

Float_t pairDZ(const muon_pog::Muon & mu1,
	       const muon_pog::Muon & mu2)
{

  // just a sanity cut as value for non filled objects is - 1000
  if ((mu1.dzbs < -990 || mu2.dzbs < -990))
    {
      std::cout << "[pairDZ]: dZ wrt BS is not valid"
		<< std::endl;
      exit(900);
    }

  Float_t result = mu1.dzbs - mu2.dzbs;

  return result;

} 

Float_t pairCos3DAngle(const TLorentzVector & muTk1,
		       const TLorentzVector & muTk2)
{

  Float_t result = cos(muTk1.Angle(muTk2.Vect()));

  return result;

} 

Float_t pTBalance(const TLorentzVector & muTk1,
		  const TLorentzVector & muTk2)
{

  Float_t ptMax = std::max(muTk1.Pt(),muTk2.Pt());
  Float_t ptMin = std::min(muTk1.Pt(),muTk2.Pt());
  

  Float_t result = ptMax / ptMin;

  return result;

} 


#endif


