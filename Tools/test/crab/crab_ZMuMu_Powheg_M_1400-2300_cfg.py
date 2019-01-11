from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'ZMuMu_Powheg_M_1400-2300_v2'

config.section_('JobType')
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = '../muonPogNtuples_cfg.py'
#config.JobType.outputFiles = ['muonNTuple.root']
config.JobType.pyCfgParams = ['globalTag=94X_mc2017_realistic_MuonTrackFix_01',
                              'ntupleName=muonPOGNtuple.root',
                              'hasRaw=True',
                              'nEvents=-1',
                              'runOnMC=True',
                              'hltPathFilter=all',
                              'minMuPt=45.0',
                              'minNMu=2'
               ]
config.JobType.allowUndistributedCMSSW = True  # To fix cmssw releases

config.section_('Data')
config.Data.inputDataset = '/ZToMuMu_NNPDF31_13TeV-powheg_M_1400_2300/RunIIFall17DRPremix-MUOTrackFix_94X_mc2017_realistic_MuonTrackFix_01_ext1-v1/AODSIM'

config.Data.useParent = True

config.Data.splitting    = 'FileBased'
config.Data.unitsPerJob  = 10  # Since files based, 10 files per job
config.Data.inputDBS     = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
config.Data.outLFNDirBase  = '/store/user/battilan/NTuplesMuonPOG/ZMu/'

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

