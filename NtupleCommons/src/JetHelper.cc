/*
 * JetHelper.cc
 *
 *  Created on: Jan 27, 2017
 *      Author: hqu
 */

#include "DeepNTuples/NtupleCommons/interface/JetHelper.h"

namespace deepntuples {

JetHelper::JetHelper(const pat::Jet* jet, const edm::Handle<reco::CandidateView> &pfcands) : jet_(jet) {
  if (!jet) throw cms::Exception("[JetHelper::JetHelper] Null pointer for input jet!");
  if (jet->nSubjetCollections() == 0) throw cms::Exception("[JetHelper::JetHelper] No subjet collection for input jet!");
  initializeConstituents(pfcands);
}

void JetHelper::initializeConstituents(const edm::Handle<reco::CandidateView> &pfcands) {
  subjets_.clear();
  uncorr_subjets_.clear();
  daughters_.clear();
  puppi_wgt_cache_.clear();

  // get subjets
  auto subjets = jet_->subjets();
  for (const auto &sj : subjets){
    subjets_.push_back(&(*sj));
    uncorr_subjets_.push_back(&(*sj));
  }
  // sort subjets by pt
  std::sort(subjets_.begin(), subjets_.end(),
      [](const pat::Jet* p1, const pat::Jet* p2){return p1->pt()>p2->pt();});

  std::sort(uncorr_subjets_.begin(), uncorr_subjets_.end(),
      [](const pat::Jet* p1, const pat::Jet* p2){return p1->correctedP4("Uncorrected").pt()>p2->correctedP4("Uncorrected").pt();});

  // get all consitituents
  for (unsigned idau=0; idau<jet_->numberOfDaughters(); ++idau){
    auto dauPtr = jet_->daughterPtr(idau);
    if (dauPtr->numberOfDaughters()>0){
      // is a subjet
      const auto *sj = dynamic_cast<const pat::Jet*>(&(*dauPtr));
      if (!sj) throw cms::Exception("[JetHelper::initializeConstituents] Cannot convert to subjet!");
      // add all daughters
      for (unsigned k=0; k<sj->numberOfDaughters(); ++k){
        const auto& candPtr = sj->daughterPtr(k);
        const auto *cand = dynamic_cast<const pat::PackedCandidate*>(&(*candPtr));
        if (cand->puppiWeight() < 0.01) continue; // [94X] ignore particles w/ extremely low puppi weights
        // Here we get the original PackedCandidate as stored in MiniAOD (i.e., not puppi weighted)
        // https://github.com/cms-sw/cmssw/pull/28035
        daughters_.push_back(pfcands->ptrAt(dauPtr.key()));
        // For the Puppi weight, we get it from the new candidate in case it is recomputed
        puppi_wgt_cache_[dauPtr.key()] = cand->puppiWeight();
      }
    }else{
      const auto& candPtr = dauPtr;
      const auto *cand = dynamic_cast<const pat::PackedCandidate*>(&(*candPtr));
      if (cand->puppiWeight() < 0.01) continue; // [94X] ignore particles w/ extremely low puppi weights
      // Here we get the original PackedCandidate as stored in MiniAOD (i.e., not puppi weighted)
      // https://github.com/cms-sw/cmssw/pull/28035
      daughters_.push_back(pfcands->ptrAt(dauPtr.key()));
      // For the Puppi weight, we get it from the new candidate in case it is recomputed
      puppi_wgt_cache_[dauPtr.key()] = cand->puppiWeight();
    }
  }
  // sort by original pt
  std::sort(daughters_.begin(), daughters_.end(), [](const reco::CandidatePtr &a, const reco::CandidatePtr &b){ return a->pt() > b->pt(); });

}

std::pair<double, double> JetHelper::getCorrectedPuppiSoftDropMass(const std::vector<const pat::Jet*> &puppisubjets) const {
  double sdpuppimass = 0;
  if (puppisubjets.size()==1){
    sdpuppimass = puppisubjets[0]->correctedP4(0).mass();
  }else if (puppisubjets.size()>=2){
    sdpuppimass = (puppisubjets[0]->correctedP4(0) + puppisubjets[1]->correctedP4(0)).mass();
  }
  double pt = jet_->pt();
  double eta = jet_->eta();
  double gencorr = 1.006261 + ((-1.061605) * pow(pt*0.079990,-1.204538));
  double recocorr = 1;
  if (std::abs(eta) <= 1.3){
    recocorr = 1.093020+(-0.000150068)*pt+(3.44866e-07)*pow(pt,2)+(-2.68100e-10)*pow(pt,3)+(8.67440e-14)*pow(pt,4)+(-1.00114e-17)*pow(pt,5);
  }else{
    recocorr = 1.272115+(-0.000571640)*pt+(8.37289e-07)*pow(pt,2)+(-5.20433e-10)*pow(pt,3)+(1.45375e-13)*pow(pt,4)+(-1.50389e-17)*pow(pt,5);
  }
  return std::make_pair(sdpuppimass, sdpuppimass*gencorr*recocorr);
}

std::vector<math::XYZTLorentzVector> JetHelper::getLBGStransformedConstituents(double jetRescale_m0, double jetLorentzBoost_e0) const
{
    /*
    Note: do not use jet().<properties> as they can be different from the sum of the jet constituents, due to weights (such as puppi).
    */
    
    double jetLorentzBoost_p0 = std::sqrt(jetLorentzBoost_e0*jetLorentzBoost_e0 - jetRescale_m0*jetRescale_m0);
    
    fastjet::Strategy fj_strategy = fastjet::Best;
    fastjet::RecombinationScheme fj_recombScheme = fastjet::E_scheme;
    
    fastjet::JetDefinition fj_fatJetExcSubJetDef(fastjet::kt_algorithm, 1.0, fj_recombScheme, fj_strategy);
    fastjet::JetDefinition fj_akJetReclusterDef(fastjet::antikt_algorithm, 1000.0, fj_recombScheme, fj_strategy);
    
    const auto& jetConstis = getJetConstituents();
    int nConsti = jetConstis.size();
    
    std::vector<fastjet::PseudoJet> v_fjInput;
    
    for(int iConsti = 0; iConsti < nConsti; iConsti++)
    {
        const auto* packed_cand = dynamic_cast<const pat::PackedCandidate*>(jetConstis.at(iConsti).get());
        
        v_fjInput.push_back(fastjet::PseudoJet(
            packed_cand->px(),
            packed_cand->py(),
            packed_cand->pz(),
            packed_cand->energy()
        ));
    }
    
    fastjet::ClusterSequence fj_jet_clustSeq(v_fjInput, fj_akJetReclusterDef);
    std::vector <fastjet::PseudoJet> fj_jets = fj_jet_clustSeq.inclusive_jets();
    fastjet::PseudoJet fj_reclusJet = fj_jets.at(0);
    
    //printf(
    //    "jet: e %f, m %f, "
    //    "\n",
    //    jet().energy(),
    //    jet().mass()
    //);
    
    //printf(
    //    "fj_reclusJet: e %f, m %f, "
    //    "\n",
    //    fj_reclusJet.e(),
    //    fj_reclusJet.m()
    //);
    
    fastjet::ClusterSequence fj_jetExcSubJet_clustSeq(fj_reclusJet.constituents(), fj_fatJetExcSubJetDef);
    std::vector <fastjet::PseudoJet> fj_jetExcSubJets = fj_jetExcSubJet_clustSeq.exclusive_jets_up_to(3);
    fj_jetExcSubJets = sorted_by_E(fj_jetExcSubJets);
    int nFatJetExcSubJet = fj_jetExcSubJets.size();
    
    std::vector <math::XYZTLorentzVector> v_jetTranformAxis;
    
    v_jetTranformAxis.push_back(math::XYZTLorentzVector(
        fj_reclusJet.px(),
        fj_reclusJet.py(),
        fj_reclusJet.pz(),
        fj_reclusJet.e()
    ));
    
    if(nFatJetExcSubJet >= 2)
    {
        v_jetTranformAxis.push_back(math::XYZTLorentzVector(
            fj_jetExcSubJets.at(0).px(),
            fj_jetExcSubJets.at(0).py(),
            fj_jetExcSubJets.at(0).pz(),
            fj_jetExcSubJets.at(0).e()
        ));
    }
    
    else
    {
        v_jetTranformAxis.push_back(math::XYZTLorentzVector(0, 0, 0, 0));
    }
    
    
    if(nFatJetExcSubJet >= 3)
    {
        v_jetTranformAxis.push_back(math::XYZTLorentzVector(
            fj_jetExcSubJets.at(1).px(),
            fj_jetExcSubJets.at(1).py(),
            fj_jetExcSubJets.at(1).pz(),
            fj_jetExcSubJets.at(1).e()
        ));
    }
    
    else
    {
        v_jetTranformAxis.push_back(math::XYZTLorentzVector(0, 0, 0, 0));
    }
    
    
    // Find the GS axes
    std::vector <math::XYZVector> v_GSaxis;
    
    for(int iAxis = 0; iAxis < (int) v_jetTranformAxis.size(); iAxis++)
    {
        math::XYZVector axis = v_jetTranformAxis.at(iAxis).Vect().Unit();
        
        for(int iGSaxis = 0; iGSaxis < (int) v_GSaxis.size(); iGSaxis++)
        {
            axis -= axis.Dot(v_GSaxis.at(iGSaxis)) * v_GSaxis.at(iGSaxis);
        }
        
        axis = axis.unit();
        v_GSaxis.push_back(axis);
    }
    
    math::XYZTLorentzVector jet_4mom_tr(0, 0, 0, 0);
    std::vector <math::XYZTLorentzVector> v_consti_4mom_tr(nConsti);
    
    // Transform the 3 momenta
    for(int iConsti = 0; iConsti < nConsti; iConsti++)
    {
        const auto* packed_cand = dynamic_cast<const pat::PackedCandidate*>(jetConstis.at(iConsti).get());
        math::XYZTLorentzVector consti_4mom = packed_cand->p4();
        
        std::vector<double> consti_3mom_tr(v_GSaxis.size());
        math::XYZTLorentzVector consti_4mom_tr;
        
        for(int iGSaxis = 0; iGSaxis < (int) v_GSaxis.size(); iGSaxis++)
        {
            consti_3mom_tr.at(iGSaxis) = consti_4mom.Vect().Dot(v_GSaxis.at(iGSaxis));
        }
        
        consti_4mom_tr.SetPxPyPzE(consti_3mom_tr.at(0), consti_3mom_tr.at(1), consti_3mom_tr.at(2), consti_4mom.e());
        v_consti_4mom_tr.at(iConsti) = consti_4mom_tr;
        jet_4mom_tr += consti_4mom_tr;
    }
    
    // Rescale
    // Very rarely, the mass can be a very small -ve value. Hence use fabs().
    double rescaleFactor = jetRescale_m0 / std::fabs(jet_4mom_tr.mass());
    jet_4mom_tr *= rescaleFactor;
    
    double boostDir = -1;
    
    if(jet_4mom_tr.e() < jetLorentzBoost_e0)
    {
        boostDir *= -1;
    }
    
    // Boost
    double boostGamma = 1.0/(jetRescale_m0*jetRescale_m0) * (
        jet_4mom_tr.e() * jetLorentzBoost_e0 -
        jetLorentzBoost_p0 * std::sqrt(jet_4mom_tr.Vect().Mag2())
    );
    
    double boostBeta = boostDir * std::sqrt(1.0 - 1.0/(boostGamma*boostGamma));
    jet_4mom_tr.SetPxPyPzE(0, 0, 0, 0);
    
    // Rescale and boost the constituents
    for(int iConsti = 0; iConsti < nConsti; iConsti++)
    {
        v_consti_4mom_tr.at(iConsti) *= rescaleFactor;
        v_consti_4mom_tr.at(iConsti) = ROOT::Math::VectorUtil::boostX(v_consti_4mom_tr.at(iConsti), boostBeta);
        
        jet_4mom_tr += v_consti_4mom_tr.at(iConsti);
    }
    
    return v_consti_4mom_tr;
}

} /* namespace deepntuples */