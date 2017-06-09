from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'SingleMuonRun2016H_PromptReco_v2'

config.section_('JobType')
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = '../muonPogNtuples_cfg.py'
#config.JobType.outputFiles = ['muonNTuple.root']
config.JobType.pyCfgParams = ['globalTag=80X_dataRun2_Prompt_v9',
                              'ntupleName=muonPOGNtuple_SingleMuonRun2016H_v2_PromptReco.root',
                              'nEvents=-1',
                              'runOnMC=False',
                              'hltPathFilter=all',
                              'minMuPt=45.0',
                              'minNMu=2'
               ]
config.JobType.allowUndistributedCMSSW = True  # To fix cmssw releases

config.section_('Data')
config.Data.inputDataset = '/SingleMuon/Run2016H-PromptReco-v2/AOD'

config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'

config.Data.splitting    = 'LumiBased'
config.Data.unitsPerJob  = 150  # Since files based, 10 files per job
config.Data.inputDBS     = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
config.Data.outLFNDirBase  = '/store/user/battilan/NTuplesMuonPOG/'

config.section_('Site')
config.Site.storageSite = 'T3_IT_Bologna'

