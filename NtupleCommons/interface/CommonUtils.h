#ifndef NTUPLECOMMONS_INTERFACE_COMMONUTILS_H_
#define NTUPLECOMMONS_INTERFACE_COMMONUTILS_H_

#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Vector.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include <algorithm>
#include <map>
#include <unordered_set>

#include "TString.h"
#include "Math/VectorUtil.h"

namespace commonutils {
    float get_cosThetaStar(
        math::XYZTLorentzVector p4_top,
        math::XYZTLorentzVector p4_b,
        math::XYZTLorentzVector p4_dau1,
        math::XYZTLorentzVector p4_dau2,
        int topType // 0: hadronic, 1: leptonic
    );
}

#endif