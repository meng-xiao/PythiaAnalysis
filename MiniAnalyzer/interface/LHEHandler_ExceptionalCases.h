// Special PDF cases

// madgraph_0offset special case with NNPDF31_lo_as_0130 as the PDF choice
bool test_specialPDF_NNPDF31_LO_as_0130_Madgraph_0offset_Case1();
// madgraph_1000offset special case with POWHEG-like numbering and, unfortunately, NNPDF30_nlo_nf_4_pdfas as the PDF choice
bool test_specialPDF_NNPDF30_NLO_nf_4_pdfas_Madgraph_1000offset_POWHEGStyle_Case1();
// madgraph_1000offset special case with some numbering differences and MC-type PDF variations
bool test_specialPDF_NNPDF31_NNLO_as_0118_nf_4_Madgraph_1000offset_Case1();
// madgraph_1000offset special case with 45 muR, muF variations with 5 completely different central value choices
bool test_specialPDF_NNPDF31_NNLO_as_0118_Madgraph_1000offset_Case1();
// madgraph_1000offset special case with a single offset of PDFs at id==1011
bool test_specialPDF_NNPDF31_NNLO_as_0118_Madgraph_1000offset_Case2();
// madgraph_1000offset special case with the central PDF being 'NNPDF31_nnlo_as_0118_mc_hessian_pdfas'
bool test_specialPDF_NNPDF31_NNLO_as_0118_mc_hessian_pdfas_Madgraph_1000offset_Case1();
// madgraph_1000offset special case with some numbering differences and MC-type PDF variations
bool test_specialPDF_NNPDF31_NNLO_as_0118_nf_4_POWHEG_MadSpin_Case1();
