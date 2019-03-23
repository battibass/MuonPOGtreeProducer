from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'Giovanni_MuonFlatGun_v5_correct'

config.section_('JobType')
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = '../muonPogNtuples_cfg.py'
#config.JobType.outputFiles = ['muonNTuple.root']
config.JobType.pyCfgParams = ['globalTag=94X_mc2017_realistic_MuonTrackFix_01',
                              'ntupleName=muonPOGNtuple_FlatMuGun_v2.root',
                              'hasRaw=True',
                              'nEvents=-1',
                              'runOnMC=True',
                              'hltPathFilter=all',
                              'minMuPt=5.0',
                              'minNMu=1'
               ]
config.JobType.allowUndistributedCMSSW = True  # To fix cmssw releases

config.section_('Data')
config.Data.inputDataset = '/MuonGun_PTOT-5-2500/abbiendi-crab_MuonGun_step3_asympt_FixedGT-TrkAli2017-v3-0de73db65e1e00fe9c9194c127b87ab0/USER'

config.Data.useParent = True

config.Data.splitting    = 'FileBased'
config.Data.unitsPerJob  = 10  # Since files based, 10 files per job
config.Data.inputDBS     = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSReader/'
config.Data.outLFNDirBase  = '/store/user/battilan/NTuplesMuonPOG/FlatPGun/'

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

