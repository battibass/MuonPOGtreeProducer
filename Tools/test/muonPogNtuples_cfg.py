import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import subprocess
import sys

options = VarParsing.VarParsing()

options.register('globalTag',
                 '80X_mcRun2_asymptotic_2016_TrancheIV_v4', #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Global Tag")

options.register('nEvents',
                 -1, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Maximum number of processed events")

options.register('eosInputFolder',
                 '/store/relval/CMSSW_8_0_20/RelValZMM_13/GEN-SIM-RECO/PU25ns_80X_mcRun2_asymptotic_2016_TrancheIV_v4_Tr4GT_v4-v1/00000', #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "EOS folder with input files")

options.register('ntupleName',
                 './muonPOGNtuple_8_0_3_RelValZMM_13.root', #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Folder and name ame for output ntuple")

options.register('runOnMC',
                 True, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Run on DATA or MC")

options.register('hasRaw',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Enables DT twin mux unpacking if RAW format available")

options.register('hasMuonTagger',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Gets bad muon tagger information")

options.register('hltPathFilter',
                 "all", #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Filter on paths (now only accepts all or IsoMu20)")

options.register('minMuPt',
                 5., #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 "Skim the ntuple saving only STA || TRK || GLB muons with pT > of this value")

options.register('minNMu',
                 1, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "number of TRK or GLB muons with pT > minMuPt to pass the skim")

options.parseArguments()

if options.hltPathFilter == "all" :
    pathCut   = "all"
    filterCut = "all"
elif options.hltPathFilter == "IsoMu20" :
    pathCut   = "HLT_IsoMu20_v"
    if options.runOnMC :
        filterCut = "hltL3crIsoL1sMu16L1f0L2f10QL3f20QL3trkIsoFiltered0p09"
    else :
        filterCut = "hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p09"
        
else :
    print "[" + sys.argv[0] + "]:", "hltPathFilter=", options.hltPathFilter, "is not a valid parameter!"
    sys.exit(100)

process = cms.Process("NTUPLES")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.nEvents))

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.MagneticField_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.GlobalTag.globaltag = cms.string(options.globalTag)

process.source = cms.Source("PoolSource",
                            
        fileNames = cms.untracked.vstring(),
        secondaryFileNames = cms.untracked.vstring()

)

#files = subprocess.check_output([ "/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select", "ls", options.eosInputFolder ])
process.source.fileNames = ['file:///afs/cern.ch/user/f/federica/public/WantedMuons/Events_RunB_1.root',
                            'file:///afs/cern.ch/user/f/federica/public/WantedMuons/Events_RunB_2.root',
                            'file:///afs/cern.ch/user/f/federica/public/WantedMuons/Events_RunC.root',
                            'file:///afs/cern.ch/user/f/federica/public/WantedMuons/Events_RunD.root',
                            'file:///afs/cern.ch/user/f/federica/public/WantedMuons/Events_RunE_1.root',
                            'file:///afs/cern.ch/user/f/federica/public/WantedMuons/Events_RunE_2.root',
                            'file:///afs/cern.ch/user/f/federica/public/WantedMuons/Events_RunF.root',
                            'file:///afs/cern.ch/user/f/federica/public/WantedMuons/Events_RunG_1.root',
                            'file:///afs/cern.ch/user/f/federica/public/WantedMuons/Events_RunG_2.root',
                            'file:///afs/cern.ch/user/f/federica/public/WantedMuons/Events_RunG_3.root',
                            'file:///afs/cern.ch/user/f/federica/public/WantedMuons/Events_RunH.root'] #options.eosInputFolder+"/"+f for f in files.split() ]  

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
#process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")

from MuonPOGtreeProducer.Tools.MuonPogNtuples_cff import appendMuonPogNtuple, customiseHlt, customiseMuonCuts
    
appendMuonPogNtuple(process,options.runOnMC,"HLT",options.ntupleName,options.hasRaw,options.hasMuonTagger)

customiseHlt(process,pathCut,filterCut)
customiseMuonCuts(process,options.minMuPt,options.minNMu)


