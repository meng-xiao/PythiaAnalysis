#ifndef DaughterDataEmbedder_h
#define DaughterDataEmbedder_h

/** \class DaughterDataEmbedder
 *
 *  A set of helpers to handle userFloats and UserData
 *
 *  \author N. Amapane - CERN
 */


#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/Photon.h>
#include <DataFormats/PatCandidates/interface/PFParticle.h>
#include <vector>



#include <string>
#include <vector>
#include <sstream>


typedef pat::PFParticle        Photon;
typedef edm::Ptr<Photon>       PhotonPtr;
typedef std::vector<PhotonPtr> PhotonPtrVector;

using namespace std;

namespace userdatahelpers {

  /// Retrieve the userFloat "name" from a reco::Candidate c
  float getUserFloat(const reco::Candidate* c, const char* name);

  /// Test if the userFloat "name" from a reco::Candidate c exists
  int hasUserFloat(const reco::Candidate* c, const char* name);

  /// Add the userFloats of daughters into cand, with proper prefix 
  void embedDaughterData(pat::CompositeCandidate& cand);

  /// Add the userFloats of daughter d into cand, with proper prefix
  template <typename T>
  void embedDaughterData(pat::CompositeCandidate& cand, unsigned i, const T* d) {
    using namespace std;

    string base;
    stringstream str;
    str << "d" << i << ".";
    str >> base;
    const std::vector<string> & userLabels = d->userFloatNames();
    for (vector<string>::const_iterator name = userLabels.begin(); name!= userLabels.end(); ++name){      
      string newname = base + *name;
      cand.addUserFloat(newname, d->userFloat(*name));
    }
  }


  /// Get the daughters daughters sorted by Z1/Z2 and charge so that the order is Z1Lp, Z1Ln, Z2Lp, Z2Ln
  /// Keeps the original sorting for the same-sign collections used for CRs
  /// Use with is4l=false for a consistent behaviour with the Z+l collection
  /// Results: 
  /// labels are the prefixes for cand.userFloat() string of each lepton.
  /// fsrIndex are indices of the leptons corresponding to each photon after sorting
  void getSortedLeptons(const pat::CompositeCandidate& cand, 
			std::vector<const reco::Candidate*>& leptons, 
			std::vector<std::string>& labels, 
			std::vector<const reco::Candidate*>& fsr,
			std::vector<short>& fsrIndex,
			bool is4l=true);

  void getSortedZLeptons(const pat::CompositeCandidate& cand, 
			std::vector<const reco::Candidate*>& leptons, 
			std::vector<std::string>& labels, 
			std::vector<const reco::Candidate*>& fsr,
			std::vector<short>& fsrIndex);


}
#endif
