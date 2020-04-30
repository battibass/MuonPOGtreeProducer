# MuonPOGtreeProducer
Muon POG tree producer

## Installation instructions

```bash
cmsrel CMSSW_11_0_2 
cd CMSSW_11_0_2/src/

git clone git@github.com:battibass/MuonPOGtreeProducer.git -b 110X_powermatches

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
# For example
cmsRun muonPogNtuples_cfg.py nEvents=1000 \
       eosInputFolder=/store/group/phys_muon/SingleMuPlusPt20to2500_NoPU/step3_CMSSW_11_0_X/200403_154053//0000/ \
       ntupleName=./muonPOGNtuple_11_0_1_SingleMuMinusPt20to2500_NoPU_1k.root

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
