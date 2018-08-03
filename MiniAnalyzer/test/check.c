void check(){
	TChain *t=new TChain ("ZZTree/candTree");
	char* pt = "PROD_samples_2017_MC_ae6f8c7/HJJ0PM_M125_Chunk*/*root";
	//char* pt = "PROD_samples_2017_MC_ae6f8c7/ggH125_minloHJJ_Chunk*/*root";
	gROOT->ProcessLine(Form(".! ls %s>info.txt",pt));
	TFileCollection fc("dum","","info.txt"); 
	fc.Print();
	t->AddFileInfoList((TCollection*)fc.GetList());
	//	t->Draw("ZZMass","ZZMass>105&&ZZMass<140");
	//t->Draw("p_JJQCD_InitialQQ_SIG_ghg2_1_JHUGen_JECNominal/(p_JJQCD_InitialQQ_SIG_ghg2_1_JHUGen_JECNominal+p_JJQCD_InitialQQ_SIG_ghg4_1_JHUGen_JECNominal)","ZZMass>105&&ZZMass<140&&nCleanedJetsPt30>1");
	t->Draw("p_Gen_HJJ_SIG_ghg2_1_JHUGen/(p_Gen_HJJ_SIG_ghg2_1_JHUGen+p_Gen_HJJ_SIG_ghg4_1_JHUGen)","");
}
