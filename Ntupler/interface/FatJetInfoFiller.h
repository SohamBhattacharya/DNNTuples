/*
 * FatJetInfoFiller.h
 *
 *  Created on: May 24, 2017
 *      Author: hqu
 */

#ifndef NTUPLER_INTERFACE_FATJETINFOFILLER_H_
#define NTUPLER_INTERFACE_FATJETINFOFILLER_H_

#include "DataFormats/BTauReco/interface/ShallowTagInfo.h"
#include "DataFormats/BTauReco/interface/BoostedDoubleSVTagInfo.h"

#include "DeepNTuples/NtupleCommons/interface/CommonUtils.h"
#include "DeepNTuples/NtupleCommons/interface/NtupleBase.h"
#include "DeepNTuples/FatJetHelpers/interface/FatJetMatching.h"

#include "Math/VectorUtil.h"

namespace deepntuples {

class FatJetInfoFiller: public NtupleBase {
public:
  FatJetInfoFiller() : FatJetInfoFiller("") {}
  FatJetInfoFiller(std::string branchName, double jetR=0.8) : NtupleBase(branchName, jetR), fjmatch_(jetR, true) {}
  virtual ~FatJetInfoFiller() {}

  // get input parameters from the cfg file
  virtual void readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector &&cc) override;

  // read event content or event setup for each event
  virtual void readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;

protected:
  // declare the data branches (name, type, default values)
  virtual void book() override;
  // fill the branches
  virtual bool fill(const pat::Jet &jet, size_t jetidx, const JetHelper &jet_helper) override;

private:
  FatJetMatching fjmatch_;
  bool useReclusteredJets_ = false;
  bool isQCDSample_ = false;
  bool isTrainSample_ = false;

  bool sample_use_pythia_ = false;
  bool sample_use_herwig_ = false;
  bool sample_use_madgraph_ = false;
  bool sample_isTopLH_ = false;
  bool sample_isTopRH_ = false;

  std::string fjTagInfoName;
  std::string fjName;

  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  edm::Handle<reco::GenParticleCollection> genParticlesHandle;

};

} /* namespace deepntuples */

#endif /* NTUPLER_INTERFACE_FATJETINFOFILLER_H_ */
