/*
 * PFCompleteFiller.h
 *
 *  Created on: Sep 25, 2017
 *      Author: hqu
 */

#ifndef NTUPLER_INTERFACE_PFCOMPLETEFILLER_H_
#define NTUPLER_INTERFACE_PFCOMPLETEFILLER_H_

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DeepNTuples/BTagHelpers/interface/TrackInfoBuilder.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DeepNTuples/NtupleCommons/interface/NtupleBase.h"

#include "DeepNTuples/NtupleCommons/interface/MVAVariableHelper.h"
#include "DeepNTuples/NtupleCommons/interface/MVAVariableManager.h"

namespace deepntuples {

  class PFCompleteFiller : public NtupleBase {
  public:
    PFCompleteFiller() : PFCompleteFiller("", 0.8) {}
    PFCompleteFiller(std::string branchName, double jetR = 0.8) : NtupleBase(branchName, jetR) {}
    virtual ~PFCompleteFiller() {}

    // get input parameters from the cfg file
    virtual void readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector &&cc) override;

    // read event content or event setup for each event
    virtual void readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;

  protected:
    // declare the data branches (name, type, default values)
    virtual void book() override;
    // fill the branches
    virtual bool fill(const pat::Jet& jet, size_t jetidx, const JetHelper& jet_helper) override;

  private:
    edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transientTrackBuilderToken_;

    edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
    edm::Handle<reco::VertexCollection> vertices;

    edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> svToken_;
    edm::Handle<reco::VertexCompositePtrCandidateCollection> SVs;

    edm::ESHandle<TransientTrackBuilder> builder_;
    
    std::unique_ptr<MVAVariableHelper> mvaVarHelper_;
    
    edm::EDGetTokenT<std::vector <pat::Electron>> elToken_;
    edm::Handle<std::vector <pat::Electron>> electronsH;
    
    edm::EDGetTokenT<std::vector <pat::Muon>> muToken_;
    edm::Handle<std::vector <pat::Muon>> muonsH;
    
    std::unique_ptr<MVAVariableManager <pat::Electron>>  elMvaVarManager_;
    std::unique_ptr<MVAVariableManager <pat::Muon>>      muMvaVarManager_;
    
    std::unordered_map<std::string, std::string> m_elVarBranchName_;
    std::unordered_map<std::string, std::string> m_muVarBranchName_;
    
    std::vector <float> extraVariables_;
  };

} /* namespace deepntuples */

#endif /* NTUPLER_INTERFACE_PFCOMPLETEFILLER_H_ */
