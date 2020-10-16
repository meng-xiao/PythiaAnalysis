#include <FWCore/Utilities/interface/EDMException.h>
#include "../interface/LHEHandler.h"
// No need for LHEHandler_ExceptionalCases.h


using namespace std;


bool LHEHandler::test_specialPDF_NNPDF31_NNLO_as_0118_nf_4_Madgraph_1000offset_Case1(){
  if (year == 2016) return false;
  else if (!(year == 2017 || year == 2018)) throw cms::Exception("SpecialPDF") << "Unknown year " << year;

  for (auto const& line:LHEHeader){
    if (line.find("<weight id=\"1010\"> PDF=  320900 NNPDF31_nnlo_as_0118_nf_4 </weight>")!=std::string::npos) return true;
  }
  return false;
}
// This is a fun one...
bool LHEHandler::test_specialPDF_NNPDF31_NNLO_as_0118_Madgraph_1000offset_Case1(){
  if (year == 2016) return false;
  else if (!(year == 2017 || year == 2018)) throw cms::Exception("SpecialPDF") << "Unknown year " << year;

  // This block of strings is only part of what needs to be tested
  const std::vector<std::string> searchBlock{
    "weight id=\"1002\" MUR=\"1.0\" MUF=\"1.0\" DYN_SCALE=\"1\" PDF=\"306000\"",
    "weight id=\"1003\" MUR=\"1.0\" MUF=\"1.0\" DYN_SCALE=\"2\" PDF=\"306000\"",
    "weight id=\"1004\" MUR=\"1.0\" MUF=\"1.0\" DYN_SCALE=\"3\" PDF=\"306000\"",
    "weight id=\"1005\" MUR=\"1.0\" MUF=\"1.0\" DYN_SCALE=\"4\" PDF=\"306000\"",
    "weight id=\"1006\" MUR=\"2.0\" MUF=\"1.0\" PDF=\"306000\"",
    "weight id=\"1007\" MUR=\"2.0\" MUF=\"1.0\" DYN_SCALE=\"1\" PDF=\"306000\"",
    "weight id=\"1008\" MUR=\"2.0\" MUF=\"1.0\" DYN_SCALE=\"2\" PDF=\"306000\"",
    "weight id=\"1009\" MUR=\"2.0\" MUF=\"1.0\" DYN_SCALE=\"3\" PDF=\"306000\"",
    "weight id=\"1010\" MUR=\"2.0\" MUF=\"1.0\" DYN_SCALE=\"4\" PDF=\"306000\"",
    "weight id=\"1011\" MUR=\"0.5\" MUF=\"1.0\" PDF=\"306000\"",
    "weight id=\"1012\" MUR=\"0.5\" MUF=\"1.0\" DYN_SCALE=\"1\" PDF=\"306000\"",
    "weight id=\"1013\" MUR=\"0.5\" MUF=\"1.0\" DYN_SCALE=\"2\" PDF=\"306000\"",
    "weight id=\"1014\" MUR=\"0.5\" MUF=\"1.0\" DYN_SCALE=\"3\" PDF=\"306000\"",
    "weight id=\"1015\" MUR=\"0.5\" MUF=\"1.0\" DYN_SCALE=\"4\" PDF=\"306000\"",
    "weight id=\"1016\" MUR=\"1.0\" MUF=\"2.0\" PDF=\"306000\"",
    "weight id=\"1017\" MUR=\"1.0\" MUF=\"2.0\" DYN_SCALE=\"1\" PDF=\"306000\"",
    "weight id=\"1018\" MUR=\"1.0\" MUF=\"2.0\" DYN_SCALE=\"2\" PDF=\"306000\"",
    "weight id=\"1019\" MUR=\"1.0\" MUF=\"2.0\" DYN_SCALE=\"3\" PDF=\"306000\"",
    "weight id=\"1020\" MUR=\"1.0\" MUF=\"2.0\" DYN_SCALE=\"4\" PDF=\"306000\"",
    "weight id=\"1021\" MUR=\"2.0\" MUF=\"2.0\" PDF=\"306000\"",
    "weight id=\"1022\" MUR=\"2.0\" MUF=\"2.0\" DYN_SCALE=\"1\" PDF=\"306000\"",
    "weight id=\"1023\" MUR=\"2.0\" MUF=\"2.0\" DYN_SCALE=\"2\" PDF=\"306000\"",
    "weight id=\"1024\" MUR=\"2.0\" MUF=\"2.0\" DYN_SCALE=\"3\" PDF=\"306000\"",
    "weight id=\"1025\" MUR=\"2.0\" MUF=\"2.0\" DYN_SCALE=\"4\" PDF=\"306000\"",
    "weight id=\"1026\" MUR=\"0.5\" MUF=\"2.0\" PDF=\"306000\"",
    "weight id=\"1027\" MUR=\"0.5\" MUF=\"2.0\" DYN_SCALE=\"1\" PDF=\"306000\"",
    "weight id=\"1028\" MUR=\"0.5\" MUF=\"2.0\" DYN_SCALE=\"2\" PDF=\"306000\"",
    "weight id=\"1029\" MUR=\"0.5\" MUF=\"2.0\" DYN_SCALE=\"3\" PDF=\"306000\"",
    "weight id=\"1030\" MUR=\"0.5\" MUF=\"2.0\" DYN_SCALE=\"4\" PDF=\"306000\"",
    "weight id=\"1031\" MUR=\"1.0\" MUF=\"0.5\" PDF=\"306000\"",
    "weight id=\"1032\" MUR=\"1.0\" MUF=\"0.5\" DYN_SCALE=\"1\" PDF=\"306000\"",
    "weight id=\"1033\" MUR=\"1.0\" MUF=\"0.5\" DYN_SCALE=\"2\" PDF=\"306000\"",
    "weight id=\"1034\" MUR=\"1.0\" MUF=\"0.5\" DYN_SCALE=\"3\" PDF=\"306000\"",
    "weight id=\"1035\" MUR=\"1.0\" MUF=\"0.5\" DYN_SCALE=\"4\" PDF=\"306000\"",
    "weight id=\"1036\" MUR=\"2.0\" MUF=\"0.5\" PDF=\"306000\"",
    "weight id=\"1037\" MUR=\"2.0\" MUF=\"0.5\" DYN_SCALE=\"1\" PDF=\"306000\"",
    "weight id=\"1038\" MUR=\"2.0\" MUF=\"0.5\" DYN_SCALE=\"2\" PDF=\"306000\"",
    "weight id=\"1039\" MUR=\"2.0\" MUF=\"0.5\" DYN_SCALE=\"3\" PDF=\"306000\"",
    "weight id=\"1040\" MUR=\"2.0\" MUF=\"0.5\" DYN_SCALE=\"4\" PDF=\"306000\"",
    "weight id=\"1041\" MUR=\"0.5\" MUF=\"0.5\" PDF=\"306000\"",
    "weight id=\"1042\" MUR=\"0.5\" MUF=\"0.5\" DYN_SCALE=\"1\" PDF=\"306000\"",
    "weight id=\"1043\" MUR=\"0.5\" MUF=\"0.5\" DYN_SCALE=\"2\" PDF=\"306000\"",
    "weight id=\"1044\" MUR=\"0.5\" MUF=\"0.5\" DYN_SCALE=\"3\" PDF=\"306000\"",
    "weight id=\"1045\" MUR=\"0.5\" MUF=\"0.5\" DYN_SCALE=\"4\" PDF=\"306000\""
  };

  bool matchAll=true;
  bool beginTest=false;
  std::vector<std::string>::const_iterator itBlock=searchBlock.cbegin();
  for (auto const& line:LHEHeader){
    if (itBlock==searchBlock.cend()) break;
    if (!beginTest && line.find(*itBlock)!=std::string::npos){
      beginTest=true;
      itBlock++;
    }
    else if (beginTest){
      matchAll &= (line.find(*itBlock)!=std::string::npos);
      if (!matchAll) break;
      itBlock++;
    }
  }

  return (beginTest && matchAll);
}
// This is just a waste of storage...
bool LHEHandler::test_specialPDF_NNPDF31_NNLO_as_0118_Madgraph_1000offset_Case2(){
  if (year == 2016) return false;
  else if (!(year == 2017 || year == 2018)) throw cms::Exception("SpecialPDF") << "Unknown year " << year;

  const std::vector<std::vector<std::string>> searchBlock{
    { "weight id", "1001", "muR=0.10000E+01 muF=0.10000E+01" },
    { "weight id", "1011", "MUR=\"1.0\" MUF=\"1.0\"" }
  };

  std::vector<std::vector<std::string>>::const_iterator itBlock=searchBlock.cbegin();
  for (auto const& line:LHEHeader){
    if (itBlock==searchBlock.cend()) break;
    bool foundAll=true;
    for (std::string const& sline:(*itBlock)) foundAll &= (line.find(sline)!=std::string::npos);
    if (foundAll) itBlock++;
  }
  return (!searchBlock.empty() && itBlock==searchBlock.cend());
}
// This is a straightforward one...
bool LHEHandler::test_specialPDF_NNPDF31_NNLO_as_0118_mc_hessian_pdfas_Madgraph_1000offset_Case1(){
  if (year == 2016) return false;
  else if (!(year == 2017 || year == 2018)) throw cms::Exception("SpecialPDF") << "Unknown year " << year;

  const std::vector<std::vector<std::string>> searchBlock{
    { "weight", "325300", "1001" }, // NNPDF31_nnlo_as_0118_mc_hessian_pdfas central PDF

    { "weight", "325300", "1010" }, // NNPDF31_nnlo_as_0118_mc_hessian_pdfas weight vars. begin
    { "weight", "325402", "1112" }, // NNPDF31_nnlo_as_0118_mc_hessian_pdfas end

    { "weight", "316200", "1113" }, // NNPDF31_nnlo_as_0118_mc begin
    { "weight", "316300", "1213" }, // NNPDF31_nnlo_as_0118_mc end

    { "weight", "306000", "1214" }, // NNPDF31_nnlo_hessian_pdfas begin
    { "weight", "306102", "1316" }, // NNPDF31_nnlo_hessian_pdfas end

    { "weight", "305800", "1325" }, // NNPDF31_nlo_hessian_pdfas begin
    { "weight", "305902", "1427" }, // NNPDF31_nlo_hessian_pdfas end

    { "weight", "292200", "2177" }, // NNPDF30_nlo_nf_5_pdfas begin
    { "weight", "292302", "2279" }, // NNPDF30_nlo_nf_5_pdfas end

    { "weight", "292600", "2280" }, // NNPDF30_nnlo_nf_5_pdfas
    { "weight", "315000", "2281" }, // NNPDF31_lo_as_0118
    { "weight", "315200", "2282" }, // NNPDF31_lo_as_0130
    { "weight", "262000", "2283" }, // NNPDF30_lo_as_0118
    { "weight", "263000", "2284" } // NNPDF30_lo_as_0130
  };

  std::vector<std::vector<std::string>>::const_iterator itBlock=searchBlock.cbegin();
  for (auto const& line:LHEHeader){
    if (itBlock==searchBlock.cend()) break;
    bool foundAll=true;
    for (std::string const& sline:(*itBlock)) foundAll &= (line.find(sline)!=std::string::npos);
    if (foundAll) itBlock++;
  }
  return (!searchBlock.empty() && itBlock==searchBlock.cend());
}

bool LHEHandler::test_specialPDF_NNPDF30_NLO_nf_4_pdfas_Madgraph_1000offset_POWHEGStyle_Case1(){
  if (year == 2016) return false;
  else if (!(year == 2017 || year == 2018)) throw cms::Exception("SpecialPDF") << "Unknown year " << year;

  const std::vector<std::vector<std::string>> searchBlock{
    {"weight id", "1001", "muR=0.10000E+01 muF=0.10000E+01"},
    {"weight id", "2001", "pdfset=292001"}
  };

  std::vector<std::vector<std::string>>::const_iterator itBlock=searchBlock.cbegin();
  for (auto const& line:LHEHeader){
    if (itBlock==searchBlock.cend()) break;
    bool foundAll=true;
    for (std::string const& sline:(*itBlock)) foundAll &= (line.find(sline)!=std::string::npos);
    if (foundAll) itBlock++;
  }
  return (!searchBlock.empty() && itBlock==searchBlock.cend());
}
// This is a case where the ordering is sort of wrong...
bool LHEHandler::test_specialPDF_NNPDF31_NNLO_as_0118_nf_4_POWHEG_MadSpin_Case1(){
  if (year == 2016) return false;
  else if (!(year == 2017 || year == 2018)) throw cms::Exception("SpecialPDF") << "Unknown year " << year;

  const std::vector<std::vector<std::string>> searchBlock{
    { "weight id", "1001", "lhapdf=320900" },
    { "weight id", "1500", "lhapdf=306000" },
    { "weight id", "2000", "lhapdf=320900" }
  };

  std::vector<std::vector<std::string>>::const_iterator itBlock=searchBlock.cbegin();
  for (auto const& line:LHEHeader){
    if (itBlock==searchBlock.cend()) break;
    bool foundAll=true;
    for (std::string const& sline:(*itBlock)) foundAll &= (line.find(sline)!=std::string::npos);
    if (foundAll) itBlock++;
  }
  return (!searchBlock.empty() && itBlock==searchBlock.cend());
}
bool LHEHandler::test_specialPDF_NNPDF31_LO_as_0130_Madgraph_0offset_Case1(){
  if (year == 2016) return false;
  else if (!(year == 2017 || year == 2018)) throw cms::Exception("SpecialPDF") << "Unknown year " << year;

  const std::vector<std::vector<std::string>> searchBlock{
    { "id", "1", "315200" }, // NNPDF31_lo_as_0130
    { "id", "111", "306000" }, // NNPDF31_nnlo_hessian_pdfas
    { "id", "222", "305800" }, // NNPDF31_nlo_hessian_pdfas
    { "id", "1074", "292200" } // NNPDF30_nlo_nf_5_pdfas
  };

  std::vector<std::vector<std::string>>::const_iterator itBlock=searchBlock.cbegin();
  for (auto const& line:LHEHeader){
    if (itBlock==searchBlock.cend()) break;
    bool foundAll=true;
    for (std::string const& sline:(*itBlock)) foundAll &= (line.find(sline)!=std::string::npos);
    if (foundAll) itBlock++;
  }
  return (!searchBlock.empty() && itBlock==searchBlock.cend());
}
