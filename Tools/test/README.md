# Instructions to produce ntuples and DataFrame for HighPt ML studies

## Ntuple production

The main cfg to be used for ntuple production `muonPogNtuples_cfg.py` is configurable.

From this directory, you can inspect the available configuration parameters by issuing:

```bash
python muonPogNtuples_cfg.py --print # this will give you the default input parameters of the filler. 
```

In particular:
* you should be able to change the input folder with EDM files to be processed: `eosInputFolder`
  * the particle guns to be used for this study are in:
    * `/store/group/phys_muon/SingleMuPlusPt20to2500_NoPU/step3_CMSSW_11_0_X/200403_154053//0000/`
    * `/store/group/phys_muon/SingleMuMinusPt20to2500_NoPU/step3_CMSSW_11_0_X/200403_153823/0000/
* you should be able to control the number of processed events: `nEvents`
  * NOTE: `nEvents=-1` means process ALL events
* you should be able to change the output root file name: `ntupleName`
  * E.g: `muonPOGNtuple_11_0_1_SingleMuPlusPt20to2500_NoPU_50k.root`

## DataFrame production (and storage as CSV)

The main cfg to be used for the production of a flat table is `muon_data_frame_creator.py`

From this directory, you can inspect the available configuration parameters by issuing:

```bash
python muon_data_frame_creator.py --help
```

It allows to dynamically include in the DataFrame variables that are in the ntuple (e.g. gen particle pt).

Additionaly it is setup to generate a few quantities that need a "non trivial" logic to be computed (e.g. track - segment pulls).

Its main configuration parameters are:
1. `-n, --nEvents` : the number of events to be processed
1. `-g, --gen` : a string with the parameters from GEN particles to be included in the data frame (enclosed in "", as: "pt eta")
1. `-m, --mu` : a string with the parameters from RECO muons to be included in the data frame (enclosed in "", as: "pt eta")
1. `-r, --refit` : a string with the parameters from RECO muon REFITs to be included in the data frame (enclosed in "", as: "pt eta")

## Ntuple format notes:

The interface of muon Ntuples is defined in : 
1. `MuonPOG/Tools/src/MuonPogTree.h` for ntuple production
2. `MuonPOG/Tools/test/interface/MuonPogTree.h` for DataFrame creation
they should be kept in synch (e.g. after a change of the first, a change on the secon should follow.
