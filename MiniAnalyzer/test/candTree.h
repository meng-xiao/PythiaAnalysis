//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Aug  1 18:20:03 2018 by ROOT version 6.06/00
// from TTree candTree/Event Summary
// found on file: ZZ4lAnalysis.root
//////////////////////////////////////////////////////////

#ifndef candTree_h
#define candTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class candTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           RunNumber;
   Long64_t        EventNumber;
   Int_t           LumiNumber;
   Short_t         NRecoMu;
   Short_t         NRecoEle;
   Short_t         Nvtx;
   Short_t         NObsInt;
   Float_t         NTrueInt;
   Float_t         PFMET;
   Short_t         nCleanedJets;
   Short_t         nCleanedJetsPt30;
   Short_t         trigWord;
   Short_t         evtPassMETFilter;
   Float_t         ZZMass;
   Short_t         ZZsel;
   Float_t         ZZPt;
   Float_t         ZZEta;
   Float_t         ZZPhi;
   Int_t           CRflag;
   Float_t         Z1Mass;
   Float_t         Z1Pt;
   Short_t         Z1Flav;
   Float_t         Z2Mass;
   Float_t         Z2Pt;
   Short_t         Z2Flav;
   Float_t         costhetastar;
   Float_t         helphi;
   Float_t         helcosthetaZ1;
   Float_t         helcosthetaZ2;
   Float_t         phistarZ1;
   Float_t         phistarZ2;
   Float_t         xi;
   Float_t         xistar;
   vector<float>   *LepPt;
   vector<float>   *LepEta;
   vector<float>   *LepPhi;
   vector<short>   *LepLepId;
   vector<float>   *JetPt;
   vector<float>   *JetEta;
   vector<float>   *JetPhi;
   vector<float>   *JetMass;
   Float_t         DiJetMass;
   Float_t         DiJetDEta;
   Short_t         nExtraLep;
   Short_t         nExtraZ;
   vector<float>   *ExtraLepPt;
   vector<float>   *ExtraLepEta;
   vector<float>   *ExtraLepPhi;
   vector<short>   *ExtraLepLepId;
   Float_t         ZXFakeweight;
   Short_t         genFinalState;
   Int_t           genProcessId;
   Float_t         genHEPMCweight;
   Float_t         genHEPMCweight_NNLO;
   Float_t         genHEPMCweight_POWHEGonly;
   Float_t         PUWeight;
   Float_t         dataMCWeight;
   Float_t         trigEffWeight;
   Float_t         overallEventWeight;
   Float_t         HqTMCweight;
   Float_t         xsec;
   Float_t         genxsec;
   Float_t         genBR;
   Short_t         genExtInfo;
   Float_t         GenHMass;
   Float_t         GenHPt;
   Float_t         GenHRapidity;
   Float_t         GenZ1Mass;
   Float_t         GenZ1Pt;
   Float_t         GenZ1Phi;
   Float_t         GenZ1Flav;
   Float_t         GenZ2Mass;
   Float_t         GenZ2Pt;
   Float_t         GenZ2Phi;
   Float_t         GenZ2Flav;
   Float_t         GenLep1Pt;
   Float_t         GenLep1Eta;
   Float_t         GenLep1Phi;
   Short_t         GenLep1Id;
   Float_t         GenLep2Pt;
   Float_t         GenLep2Eta;
   Float_t         GenLep2Phi;
   Short_t         GenLep2Id;
   Float_t         GenLep3Pt;
   Float_t         GenLep3Eta;
   Float_t         GenLep3Phi;
   Short_t         GenLep3Id;
   Float_t         GenLep4Pt;
   Float_t         GenLep4Eta;
   Float_t         GenLep4Phi;
   Short_t         GenLep4Id;
   Float_t         GenAssocLep1Pt;
   Float_t         GenAssocLep1Eta;
   Float_t         GenAssocLep1Phi;
   Short_t         GenAssocLep1Id;
   Float_t         GenAssocLep2Pt;
   Float_t         GenAssocLep2Eta;
   Float_t         GenAssocLep2Phi;
   Short_t         GenAssocLep2Id;
   vector<float>   *LHEMotherPz;
   vector<float>   *LHEMotherE;
   vector<short>   *LHEMotherId;
   vector<float>   *LHEDaughterPt;
   vector<float>   *LHEDaughterEta;
   vector<float>   *LHEDaughterPhi;
   vector<float>   *LHEDaughterMass;
   vector<short>   *LHEDaughterId;
   vector<float>   *LHEAssociatedParticlePt;
   vector<float>   *LHEAssociatedParticleEta;
   vector<float>   *LHEAssociatedParticlePhi;
   vector<float>   *LHEAssociatedParticleMass;
   vector<short>   *LHEAssociatedParticleId;
   Float_t         LHEPDFScale;
   Float_t         p_Gen_Dec_SIG_ghz1_1_JHUGen;
   Float_t         p_Gen_Dec_SIG_ghz4_1_JHUGen;
   Float_t         p_Gen_HJJ_SIG_ghg2_1_JHUGen;
   Float_t         p_Gen_HJJ_SIG_ghg4_1_JHUGen;

   // List of branches
   TBranch        *b_RunNumber;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_LumiNumber;   //!
   TBranch        *b_NRecoMu;   //!
   TBranch        *b_NRecoEle;   //!
   TBranch        *b_Nvtx;   //!
   TBranch        *b_NObsInt;   //!
   TBranch        *b_NTrueInt;   //!
   TBranch        *b_PFMET;   //!
   TBranch        *b_nCleanedJets;   //!
   TBranch        *b_nCleanedJetsPt30;   //!
   TBranch        *b_trigWord;   //!
   TBranch        *b_evtPassMETFilter;   //!
   TBranch        *b_ZZMass;   //!
   TBranch        *b_ZZsel;   //!
   TBranch        *b_ZZPt;   //!
   TBranch        *b_ZZEta;   //!
   TBranch        *b_ZZPhi;   //!
   TBranch        *b_CRflag;   //!
   TBranch        *b_Z1Mass;   //!
   TBranch        *b_Z1Pt;   //!
   TBranch        *b_Z1Flav;   //!
   TBranch        *b_Z2Mass;   //!
   TBranch        *b_Z2Pt;   //!
   TBranch        *b_Z2Flav;   //!
   TBranch        *b_costhetastar;   //!
   TBranch        *b_helphi;   //!
   TBranch        *b_helcosthetaZ1;   //!
   TBranch        *b_helcosthetaZ2;   //!
   TBranch        *b_phistarZ1;   //!
   TBranch        *b_phistarZ2;   //!
   TBranch        *b_xi;   //!
   TBranch        *b_xistar;   //!
   TBranch        *b_LepPt;   //!
   TBranch        *b_LepEta;   //!
   TBranch        *b_LepPhi;   //!
   TBranch        *b_LepLepId;   //!
   TBranch        *b_JetPt;   //!
   TBranch        *b_JetEta;   //!
   TBranch        *b_JetPhi;   //!
   TBranch        *b_JetMass;   //!
   TBranch        *b_DiJetMass;   //!
   TBranch        *b_DiJetDEta;   //!
   TBranch        *b_nExtraLep;   //!
   TBranch        *b_nExtraZ;   //!
   TBranch        *b_ExtraLepPt;   //!
   TBranch        *b_ExtraLepEta;   //!
   TBranch        *b_ExtraLepPhi;   //!
   TBranch        *b_ExtraLepLepId;   //!
   TBranch        *b_ZXFakeweight;   //!
   TBranch        *b_genFinalState;   //!
   TBranch        *b_genProcessId;   //!
   TBranch        *b_genHEPMCweight;   //!
   TBranch        *b_genHEPMCweight_NNLO;   //!
   TBranch        *b_genHEPMCweight_POWHEGonly;   //!
   TBranch        *b_PUWeight;   //!
   TBranch        *b_dataMCWeight;   //!
   TBranch        *b_trigEffWeight;   //!
   TBranch        *b_overallEventWeight;   //!
   TBranch        *b_HqTMCweight;   //!
   TBranch        *b_xsec;   //!
   TBranch        *b_genxsec;   //!
   TBranch        *b_genBR;   //!
   TBranch        *b_genExtInfo;   //!
   TBranch        *b_GenHMass;   //!
   TBranch        *b_GenHPt;   //!
   TBranch        *b_GenHRapidity;   //!
   TBranch        *b_GenZ1Mass;   //!
   TBranch        *b_GenZ1Pt;   //!
   TBranch        *b_GenZ1Phi;   //!
   TBranch        *b_GenZ1Flav;   //!
   TBranch        *b_GenZ2Mass;   //!
   TBranch        *b_GenZ2Pt;   //!
   TBranch        *b_GenZ2Phi;   //!
   TBranch        *b_GenZ2Flav;   //!
   TBranch        *b_GenLep1Pt;   //!
   TBranch        *b_GenLep1Eta;   //!
   TBranch        *b_GenLep1Phi;   //!
   TBranch        *b_GenLep1Id;   //!
   TBranch        *b_GenLep2Pt;   //!
   TBranch        *b_GenLep2Eta;   //!
   TBranch        *b_GenLep2Phi;   //!
   TBranch        *b_GenLep2Id;   //!
   TBranch        *b_GenLep3Pt;   //!
   TBranch        *b_GenLep3Eta;   //!
   TBranch        *b_GenLep3Phi;   //!
   TBranch        *b_GenLep3Id;   //!
   TBranch        *b_GenLep4Pt;   //!
   TBranch        *b_GenLep4Eta;   //!
   TBranch        *b_GenLep4Phi;   //!
   TBranch        *b_GenLep4Id;   //!
   TBranch        *b_GenAssocLep1Pt;   //!
   TBranch        *b_GenAssocLep1Eta;   //!
   TBranch        *b_GenAssocLep1Phi;   //!
   TBranch        *b_GenAssocLep1Id;   //!
   TBranch        *b_GenAssocLep2Pt;   //!
   TBranch        *b_GenAssocLep2Eta;   //!
   TBranch        *b_GenAssocLep2Phi;   //!
   TBranch        *b_GenAssocLep2Id;   //!
   TBranch        *b_LHEMotherPz;   //!
   TBranch        *b_LHEMotherE;   //!
   TBranch        *b_LHEMotherId;   //!
   TBranch        *b_LHEDaughterPt;   //!
   TBranch        *b_LHEDaughterEta;   //!
   TBranch        *b_LHEDaughterPhi;   //!
   TBranch        *b_LHEDaughterMass;   //!
   TBranch        *b_LHEDaughterId;   //!
   TBranch        *b_LHEAssociatedParticlePt;   //!
   TBranch        *b_LHEAssociatedParticleEta;   //!
   TBranch        *b_LHEAssociatedParticlePhi;   //!
   TBranch        *b_LHEAssociatedParticleMass;   //!
   TBranch        *b_LHEAssociatedParticleId;   //!
   TBranch        *b_LHEPDFScale;   //!
   TBranch        *b_p_Gen_Dec_SIG_ghz1_1_JHUGen;   //!
   TBranch        *b_p_Gen_Dec_SIG_ghz4_1_JHUGen;   //!
   TBranch        *b_p_Gen_HJJ_SIG_ghg2_1_JHUGen;   //!
   TBranch        *b_p_Gen_HJJ_SIG_ghg4_1_JHUGen;   //!

   candTree(TTree *tree=0);
   virtual ~candTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef candTree_cxx
candTree::candTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ZZ4lAnalysis.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ZZ4lAnalysis.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("ZZ4lAnalysis.root:/ZZTree");
      dir->GetObject("candTree",tree);

   }
   Init(tree);
}

candTree::~candTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t candTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t candTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void candTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   LepPt = 0;
   LepEta = 0;
   LepPhi = 0;
   LepLepId = 0;
   JetPt = 0;
   JetEta = 0;
   JetPhi = 0;
   JetMass = 0;
   ExtraLepPt = 0;
   ExtraLepEta = 0;
   ExtraLepPhi = 0;
   ExtraLepLepId = 0;
   LHEMotherPz = 0;
   LHEMotherE = 0;
   LHEMotherId = 0;
   LHEDaughterPt = 0;
   LHEDaughterEta = 0;
   LHEDaughterPhi = 0;
   LHEDaughterMass = 0;
   LHEDaughterId = 0;
   LHEAssociatedParticlePt = 0;
   LHEAssociatedParticleEta = 0;
   LHEAssociatedParticlePhi = 0;
   LHEAssociatedParticleMass = 0;
   LHEAssociatedParticleId = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("LumiNumber", &LumiNumber, &b_LumiNumber);
   fChain->SetBranchAddress("NRecoMu", &NRecoMu, &b_NRecoMu);
   fChain->SetBranchAddress("NRecoEle", &NRecoEle, &b_NRecoEle);
   fChain->SetBranchAddress("Nvtx", &Nvtx, &b_Nvtx);
   fChain->SetBranchAddress("NObsInt", &NObsInt, &b_NObsInt);
   fChain->SetBranchAddress("NTrueInt", &NTrueInt, &b_NTrueInt);
   fChain->SetBranchAddress("PFMET", &PFMET, &b_PFMET);
   fChain->SetBranchAddress("nCleanedJets", &nCleanedJets, &b_nCleanedJets);
   fChain->SetBranchAddress("nCleanedJetsPt30", &nCleanedJetsPt30, &b_nCleanedJetsPt30);
   fChain->SetBranchAddress("trigWord", &trigWord, &b_trigWord);
   fChain->SetBranchAddress("evtPassMETFilter", &evtPassMETFilter, &b_evtPassMETFilter);
   fChain->SetBranchAddress("ZZMass", &ZZMass, &b_ZZMass);
   fChain->SetBranchAddress("ZZsel", &ZZsel, &b_ZZsel);
   fChain->SetBranchAddress("ZZPt", &ZZPt, &b_ZZPt);
   fChain->SetBranchAddress("ZZEta", &ZZEta, &b_ZZEta);
   fChain->SetBranchAddress("ZZPhi", &ZZPhi, &b_ZZPhi);
   fChain->SetBranchAddress("CRflag", &CRflag, &b_CRflag);
   fChain->SetBranchAddress("Z1Mass", &Z1Mass, &b_Z1Mass);
   fChain->SetBranchAddress("Z1Pt", &Z1Pt, &b_Z1Pt);
   fChain->SetBranchAddress("Z1Flav", &Z1Flav, &b_Z1Flav);
   fChain->SetBranchAddress("Z2Mass", &Z2Mass, &b_Z2Mass);
   fChain->SetBranchAddress("Z2Pt", &Z2Pt, &b_Z2Pt);
   fChain->SetBranchAddress("Z2Flav", &Z2Flav, &b_Z2Flav);
   fChain->SetBranchAddress("costhetastar", &costhetastar, &b_costhetastar);
   fChain->SetBranchAddress("helphi", &helphi, &b_helphi);
   fChain->SetBranchAddress("helcosthetaZ1", &helcosthetaZ1, &b_helcosthetaZ1);
   fChain->SetBranchAddress("helcosthetaZ2", &helcosthetaZ2, &b_helcosthetaZ2);
   fChain->SetBranchAddress("phistarZ1", &phistarZ1, &b_phistarZ1);
   fChain->SetBranchAddress("phistarZ2", &phistarZ2, &b_phistarZ2);
   fChain->SetBranchAddress("xi", &xi, &b_xi);
   fChain->SetBranchAddress("xistar", &xistar, &b_xistar);
   fChain->SetBranchAddress("LepPt", &LepPt, &b_LepPt);
   fChain->SetBranchAddress("LepEta", &LepEta, &b_LepEta);
   fChain->SetBranchAddress("LepPhi", &LepPhi, &b_LepPhi);
   fChain->SetBranchAddress("LepLepId", &LepLepId, &b_LepLepId);
   fChain->SetBranchAddress("JetPt", &JetPt, &b_JetPt);
   fChain->SetBranchAddress("JetEta", &JetEta, &b_JetEta);
   fChain->SetBranchAddress("JetPhi", &JetPhi, &b_JetPhi);
   fChain->SetBranchAddress("JetMass", &JetMass, &b_JetMass);
   fChain->SetBranchAddress("DiJetMass", &DiJetMass, &b_DiJetMass);
   fChain->SetBranchAddress("DiJetDEta", &DiJetDEta, &b_DiJetDEta);
   fChain->SetBranchAddress("nExtraLep", &nExtraLep, &b_nExtraLep);
   fChain->SetBranchAddress("nExtraZ", &nExtraZ, &b_nExtraZ);
   fChain->SetBranchAddress("ExtraLepPt", &ExtraLepPt, &b_ExtraLepPt);
   fChain->SetBranchAddress("ExtraLepEta", &ExtraLepEta, &b_ExtraLepEta);
   fChain->SetBranchAddress("ExtraLepPhi", &ExtraLepPhi, &b_ExtraLepPhi);
   fChain->SetBranchAddress("ExtraLepLepId", &ExtraLepLepId, &b_ExtraLepLepId);
   fChain->SetBranchAddress("ZXFakeweight", &ZXFakeweight, &b_ZXFakeweight);
   fChain->SetBranchAddress("genFinalState", &genFinalState, &b_genFinalState);
   fChain->SetBranchAddress("genProcessId", &genProcessId, &b_genProcessId);
   fChain->SetBranchAddress("genHEPMCweight", &genHEPMCweight, &b_genHEPMCweight);
   fChain->SetBranchAddress("genHEPMCweight_NNLO", &genHEPMCweight_NNLO, &b_genHEPMCweight_NNLO);
   fChain->SetBranchAddress("genHEPMCweight_POWHEGonly", &genHEPMCweight_POWHEGonly, &b_genHEPMCweight_POWHEGonly);
   fChain->SetBranchAddress("PUWeight", &PUWeight, &b_PUWeight);
   fChain->SetBranchAddress("dataMCWeight", &dataMCWeight, &b_dataMCWeight);
   fChain->SetBranchAddress("trigEffWeight", &trigEffWeight, &b_trigEffWeight);
   fChain->SetBranchAddress("overallEventWeight", &overallEventWeight, &b_overallEventWeight);
   fChain->SetBranchAddress("HqTMCweight", &HqTMCweight, &b_HqTMCweight);
   fChain->SetBranchAddress("xsec", &xsec, &b_xsec);
   fChain->SetBranchAddress("genxsec", &genxsec, &b_genxsec);
   fChain->SetBranchAddress("genBR", &genBR, &b_genBR);
   fChain->SetBranchAddress("genExtInfo", &genExtInfo, &b_genExtInfo);
   fChain->SetBranchAddress("GenHMass", &GenHMass, &b_GenHMass);
   fChain->SetBranchAddress("GenHPt", &GenHPt, &b_GenHPt);
   fChain->SetBranchAddress("GenHRapidity", &GenHRapidity, &b_GenHRapidity);
   fChain->SetBranchAddress("GenZ1Mass", &GenZ1Mass, &b_GenZ1Mass);
   fChain->SetBranchAddress("GenZ1Pt", &GenZ1Pt, &b_GenZ1Pt);
   fChain->SetBranchAddress("GenZ1Phi", &GenZ1Phi, &b_GenZ1Phi);
   fChain->SetBranchAddress("GenZ1Flav", &GenZ1Flav, &b_GenZ1Flav);
   fChain->SetBranchAddress("GenZ2Mass", &GenZ2Mass, &b_GenZ2Mass);
   fChain->SetBranchAddress("GenZ2Pt", &GenZ2Pt, &b_GenZ2Pt);
   fChain->SetBranchAddress("GenZ2Phi", &GenZ2Phi, &b_GenZ2Phi);
   fChain->SetBranchAddress("GenZ2Flav", &GenZ2Flav, &b_GenZ2Flav);
   fChain->SetBranchAddress("GenLep1Pt", &GenLep1Pt, &b_GenLep1Pt);
   fChain->SetBranchAddress("GenLep1Eta", &GenLep1Eta, &b_GenLep1Eta);
   fChain->SetBranchAddress("GenLep1Phi", &GenLep1Phi, &b_GenLep1Phi);
   fChain->SetBranchAddress("GenLep1Id", &GenLep1Id, &b_GenLep1Id);
   fChain->SetBranchAddress("GenLep2Pt", &GenLep2Pt, &b_GenLep2Pt);
   fChain->SetBranchAddress("GenLep2Eta", &GenLep2Eta, &b_GenLep2Eta);
   fChain->SetBranchAddress("GenLep2Phi", &GenLep2Phi, &b_GenLep2Phi);
   fChain->SetBranchAddress("GenLep2Id", &GenLep2Id, &b_GenLep2Id);
   fChain->SetBranchAddress("GenLep3Pt", &GenLep3Pt, &b_GenLep3Pt);
   fChain->SetBranchAddress("GenLep3Eta", &GenLep3Eta, &b_GenLep3Eta);
   fChain->SetBranchAddress("GenLep3Phi", &GenLep3Phi, &b_GenLep3Phi);
   fChain->SetBranchAddress("GenLep3Id", &GenLep3Id, &b_GenLep3Id);
   fChain->SetBranchAddress("GenLep4Pt", &GenLep4Pt, &b_GenLep4Pt);
   fChain->SetBranchAddress("GenLep4Eta", &GenLep4Eta, &b_GenLep4Eta);
   fChain->SetBranchAddress("GenLep4Phi", &GenLep4Phi, &b_GenLep4Phi);
   fChain->SetBranchAddress("GenLep4Id", &GenLep4Id, &b_GenLep4Id);
   fChain->SetBranchAddress("GenAssocLep1Pt", &GenAssocLep1Pt, &b_GenAssocLep1Pt);
   fChain->SetBranchAddress("GenAssocLep1Eta", &GenAssocLep1Eta, &b_GenAssocLep1Eta);
   fChain->SetBranchAddress("GenAssocLep1Phi", &GenAssocLep1Phi, &b_GenAssocLep1Phi);
   fChain->SetBranchAddress("GenAssocLep1Id", &GenAssocLep1Id, &b_GenAssocLep1Id);
   fChain->SetBranchAddress("GenAssocLep2Pt", &GenAssocLep2Pt, &b_GenAssocLep2Pt);
   fChain->SetBranchAddress("GenAssocLep2Eta", &GenAssocLep2Eta, &b_GenAssocLep2Eta);
   fChain->SetBranchAddress("GenAssocLep2Phi", &GenAssocLep2Phi, &b_GenAssocLep2Phi);
   fChain->SetBranchAddress("GenAssocLep2Id", &GenAssocLep2Id, &b_GenAssocLep2Id);
   fChain->SetBranchAddress("LHEMotherPz", &LHEMotherPz, &b_LHEMotherPz);
   fChain->SetBranchAddress("LHEMotherE", &LHEMotherE, &b_LHEMotherE);
   fChain->SetBranchAddress("LHEMotherId", &LHEMotherId, &b_LHEMotherId);
   fChain->SetBranchAddress("LHEDaughterPt", &LHEDaughterPt, &b_LHEDaughterPt);
   fChain->SetBranchAddress("LHEDaughterEta", &LHEDaughterEta, &b_LHEDaughterEta);
   fChain->SetBranchAddress("LHEDaughterPhi", &LHEDaughterPhi, &b_LHEDaughterPhi);
   fChain->SetBranchAddress("LHEDaughterMass", &LHEDaughterMass, &b_LHEDaughterMass);
   fChain->SetBranchAddress("LHEDaughterId", &LHEDaughterId, &b_LHEDaughterId);
   fChain->SetBranchAddress("LHEAssociatedParticlePt", &LHEAssociatedParticlePt, &b_LHEAssociatedParticlePt);
   fChain->SetBranchAddress("LHEAssociatedParticleEta", &LHEAssociatedParticleEta, &b_LHEAssociatedParticleEta);
   fChain->SetBranchAddress("LHEAssociatedParticlePhi", &LHEAssociatedParticlePhi, &b_LHEAssociatedParticlePhi);
   fChain->SetBranchAddress("LHEAssociatedParticleMass", &LHEAssociatedParticleMass, &b_LHEAssociatedParticleMass);
   fChain->SetBranchAddress("LHEAssociatedParticleId", &LHEAssociatedParticleId, &b_LHEAssociatedParticleId);
   fChain->SetBranchAddress("LHEPDFScale", &LHEPDFScale, &b_LHEPDFScale);
   fChain->SetBranchAddress("p_Gen_Dec_SIG_ghz1_1_JHUGen", &p_Gen_Dec_SIG_ghz1_1_JHUGen, &b_p_Gen_Dec_SIG_ghz1_1_JHUGen);
   fChain->SetBranchAddress("p_Gen_Dec_SIG_ghz4_1_JHUGen", &p_Gen_Dec_SIG_ghz4_1_JHUGen, &b_p_Gen_Dec_SIG_ghz4_1_JHUGen);
   fChain->SetBranchAddress("p_Gen_HJJ_SIG_ghg2_1_JHUGen", &p_Gen_HJJ_SIG_ghg2_1_JHUGen, &b_p_Gen_HJJ_SIG_ghg2_1_JHUGen);
   fChain->SetBranchAddress("p_Gen_HJJ_SIG_ghg4_1_JHUGen", &p_Gen_HJJ_SIG_ghg4_1_JHUGen, &b_p_Gen_HJJ_SIG_ghg4_1_JHUGen);
   Notify();
}

Bool_t candTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void candTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t candTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef candTree_cxx
