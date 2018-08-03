import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
		# replace 'myfile.root' with the source file you want to use
		fileNames = cms.untracked.vstring(
			#        '/store/mc/RunIISummer15GS/GluGluHToZZTo4L_M125_13TeV_powheg2_minloHJ_NNLOPS_JHUgenV702_pythia8_V2/GEN-SIM/MCRUN2_71_V1-v1/110000/02995DC1-B191-E711-AB1E-001EC9AF64B6.root '
			#'/store/mc/RunIIFall17MiniAOD/ZZTo4L_13TeV_powheg_pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v2/60000/421F9A3D-A1D7-E711-8629-0CC47A13CCFC.root'
'/store/mc/RunIISummer16MiniAODv2/GluGluHToZZTo4L_M125_13TeV_powheg2_minloHJ_NNLOPS_JHUgenV702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/44EEBBEB-E0DE-E611-A786-0CC47A07F9FE.root'
)
		)

process.demo = cms.EDAnalyzer('MiniAnalyzer',
		genjet = cms.InputTag("slimmedGenJets"),
		pruned = cms.InputTag("prunedGenParticles"),
		VVDecayMode = cms.int32(0),
		VVMode = cms.int32(1),
		failedTreeLevel = cms.int32(0),
		sampleName = cms.string('HJJ0PM_M125'),
		lheProbabilities = cms.vstring('Name:SampleProductionHypothesisJHUGen Alias:<Name> Process:SelfDefine_spin0 Production:JJQCD MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0 Cluster:BestLOAssociatedVBF isGen:1 NoBranch:1',
			'Name:HJJ_SIG_ghg2_1_JHUGen Alias:SampleProductionHypothesisJHUGen Process:SelfDefine_spin0 Production:JJQCD MatrixElement:JHUGen Couplings:ghg2=1,0 Options:DivideP=SampleProductionHypothesisJHUGen Cluster:BestLOAssociatedVBF isGen:1',
			'Name:HJJ_SIG_ghg4_1_JHUGen Process:SelfDefine_spin0 Production:JJQCD MatrixElement:JHUGen Couplings:ghg4=1,0 Options:DivideP=SampleProductionHypothesisJHUGen Cluster:BestLOAssociatedVBF isGen:1',
			)
		)


process.p = cms.Path(process.demo)






