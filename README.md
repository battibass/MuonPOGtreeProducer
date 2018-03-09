# MuonPOGtreeProducer
Muon POG tree producer

## Installation instructions

```bash
cmsrel CMSSW_9_4_4 
cd CMSSW_9_4_4/src/

git clone git@github.com:battibass/MuonPOGtreeProducer.git

cmsenv

scramv1 b -j 5
```


## Ntuples

The interface of muon Ntuples is defined in : MuonPOG/Tools/src/MuonPogTree.h

The code filling ntuples is available in : MuonPOG/Tools/plugins/MuonPogTreeProducer.cc

It fills HLT, GEN level, beam spot, vertex and muon information. It works both in AOD and miniAOD (NOTE: trigger information not filled when running in miniAOD).


To create some ntuples :

```bash
cd MuonPOG/Tools/test/
python muonPogNtuples_cfg.py --print # this will give you the default input parameters of the filler. 
                                     # As the ntuple cfg is based on VarParsing you can customise the
                                     # ntuple production via command line [1] or in a crab cfg [2] 

[1] 
cmsRun muonPogNtuples_cfg.py globalTag=80X_mcRun2_asymptotic_v5 \\
  eosInputFolder=/store/relval/CMSSW_8_0_0_patch2/RelValZMM_13/GEN-SIM-RECO/PU25ns_80X_mcRun2_asymptotic_v5_refGT-v1/10000

[2] 
https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile#CRAB_configuration_parameters (find pyCfgParams)
```

## Running on CRAB

For running trhough crab you can go to: 

    MuonPOG/Tools/test/crab/

The CRAB client can be sourced using the command below after cmsenv.

    source /cvmfs/cms.cern.ch/crab3/crab.sh
  
Check if you have writing permissions in the common area if you already asked for that. 

    crab checkwrite --site=T2_CH_CERN --lfn=/store/group/phys_muon/

For running on the most recent re-reco data you can open file example: 

    crab_SingleMuon_Run2016B-23SepReReco_Run273158.py    

whit the configuration pset parameters: 

    config.JobType.pyCfgParams = ['globalTag=80X_dataRun2_2016SeptRepro_v3',
                                  'ntupleName=muonPOGNtuple_SingleMuonRun2016B_23SepReReco.root',
                                  'nEvents=-1',
                                  'runOnMC=False',
                                  'hltPathFilter=all',
                                  'minMuPt=10.0',
                                  'minNMu=2' ]


The ntuple producer gets loaded by :

```python
from MuonPOG.Tools.MuonPogNtuples_cff import appendMuonPogNtuple
appendMuonPogNtuple(process,False,"HLT","MY_NTUPLE_NAME.root")
```

Where arguments are :

1. The CMS configuration process
2. A bool to say whether you are running on MC
3. The label of the process giving HLT results
4. The name of the output ntuple
