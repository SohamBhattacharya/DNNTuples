/*
 * PFCompleteFiller.cc
 *
 *  Created on: Sep 25, 2017
 *      Author: hqu
 */

#include <unordered_map>
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DeepNTuples/Ntupler/interface/PFCompleteFiller.h"

namespace deepntuples
{

    void PFCompleteFiller::readConfig(const edm::ParameterSet &iConfig, edm::ConsumesCollector &&cc)
    {
        transientTrackBuilderToken_ =
            cc.esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"));
        vtxToken_ = cc.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
        svToken_ = cc.consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("SVs"));

        elToken_ = cc.consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"));
        muToken_ = cc.consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"));

        mvaVarHelper_.reset(new MVAVariableHelper(cc));
        elMvaVarManager_.reset(new MVAVariableManager<pat::Electron>(iConfig.getParameter<std::string>("elMvaVariablesFile"), MVAVariableHelper::indexMap()));
        muMvaVarManager_.reset(new MVAVariableManager<pat::Muon>(iConfig.getParameter<std::string>("muMvaVariablesFile"), MVAVariableHelper::indexMap()));
    }

    void PFCompleteFiller::readEvent(const edm::Event &iEvent, const edm::EventSetup &iSetup)
    {
        iEvent.getByToken(vtxToken_, vertices);
        iEvent.getByToken(svToken_, SVs);
        
        iEvent.getByToken(elToken_, electronsH);
        iEvent.getByToken(muToken_, muonsH);
        
        extraVariables_ = mvaVarHelper_->getAuxVariables(iEvent);
        
        builder_ = iSetup.getHandle(transientTrackBuilderToken_);
    }

    void PFCompleteFiller::book()
    {
        data.add<int>("n_pfcands", 0);
        data.add<float>("npfcands", 0);

        // no puppi scaled
        data.addMulti<float>("pfcand_px");
        data.addMulti<float>("pfcand_py");
        data.addMulti<float>("pfcand_pz");
        data.addMulti<float>("pfcand_energy");

        data.addMulti<float>("pfcand_LBGS_x");
        data.addMulti<float>("pfcand_LBGS_y");
        data.addMulti<float>("pfcand_LBGS_enFrac");

        data.addMulti<float>("pfcand_pt_nopuppi");
        data.addMulti<float>("pfcand_pt_log_nopuppi");
        data.addMulti<float>("pfcand_e_log_nopuppi");

        data.addMulti<float>("pfcand_phirel");
        data.addMulti<float>("pfcand_etarel");
        data.addMulti<float>("pfcand_puppiw");
        data.addMulti<float>("pfcand_abseta");
        data.addMulti<float>("pfcand_ptWrtJet");

        data.addMulti<float>("pfcand_charge");
        data.addMulti<float>("pfcand_isMu");
        data.addMulti<float>("pfcand_isEl");
        data.addMulti<float>("pfcand_isChargedHad");
        data.addMulti<float>("pfcand_isGamma");
        data.addMulti<float>("pfcand_isNeutralHad");

        // for neutral
        data.addMulti<float>("pfcand_hcalFrac");
        data.addMulti<float>("pfcand_hcalFracCalib");

        // for charged
        data.addMulti<float>("pfcand_VTX_ass");
        data.addMulti<float>("pfcand_fromPV");
        data.addMulti<float>("pfcand_lostInnerHits");
        data.addMulti<float>("pfcand_trackHighPurity");

        // impact parameters
        data.addMulti<float>("pfcand_dz");
        data.addMulti<float>("pfcand_dzsig");
        data.addMulti<float>("pfcand_dxy");
        data.addMulti<float>("pfcand_dxysig");

        // track quality
        data.addMulti<float>("pfcand_normchi2");
        data.addMulti<float>("pfcand_quality");

        // track btag info
        // data.addMulti<float>("pfcand_btagMomentum");
        // data.addMulti<float>("pfcand_btagEta");
        data.addMulti<float>("pfcand_btagEtaRel");
        data.addMulti<float>("pfcand_btagPtRel");
        // data.addMulti<float>("pfcand_btagPPar");
        // data.addMulti<float>("pfcand_btagDeltaR");
        data.addMulti<float>("pfcand_btagPtRatio");
        data.addMulti<float>("pfcand_btagPParRatio");
        data.addMulti<float>("pfcand_btagSip3dVal");
        data.addMulti<float>("pfcand_btagSip3dSig");
        data.addMulti<float>("pfcand_btagJetDistVal");
        //  data.addMulti<float>("pfcand_btagJetDistSig"); // always gives 0?
        
        data.addMulti<float>("pfcand_isElMatched");
        data.addMulti<float>("pfcand_isMuMatched");
        
        for (int iVar = 0; iVar < elMvaVarManager_->getNVars(); iVar++)
        {
            char brName[2000];
            sprintf(brName, "pfcand_%s", elMvaVarManager_->getName(iVar).c_str());
            m_elVarBranchName_[elMvaVarManager_->getName(iVar)] = std::string(brName);
            data.addMulti<float>(brName);
        }

        for (int iVar = 0; iVar < muMvaVarManager_->getNVars(); iVar++)
        {
            char brName[2000];
            sprintf(brName, "pfcand_%s", muMvaVarManager_->getName(iVar).c_str());
            m_muVarBranchName_[muMvaVarManager_->getName(iVar)] = std::string(brName);
            data.addMulti<float>(brName);
        }
    }

    bool PFCompleteFiller::fill(const pat::Jet &jet, size_t jetidx, const JetHelper &jet_helper)
    {
        const auto &electrons = *electronsH;
        const auto &muons = *muonsH;
        
        const auto &pfCands = jet_helper.getJetConstituents();

        data.fill<int>("n_pfcands", pfCands.size());
        data.fill<float>("npfcands", pfCands.size());

        float etasign = jet.eta() > 0 ? 1 : -1;

        double jetRescale_m0 = 1.0;
        double jetLorentzBoost_e0 = 2.0;

        std::vector<math::XYZTLorentzVector> v_consti_4mom_tr = jet_helper.getLBGStransformedConstituents(jetRescale_m0, jetLorentzBoost_e0);

        int iConsti = -1;

        for (const auto &cand : pfCands)
        {
            iConsti++;
            const auto *packed_cand = dynamic_cast<const pat::PackedCandidate *>(&(*cand));

            // basic kinematics, valid for both charged and neutral
            // not puppi weighted
            data.fillMulti<float>("pfcand_px", packed_cand->px());
            data.fillMulti<float>("pfcand_py", packed_cand->py());
            data.fillMulti<float>("pfcand_pz", packed_cand->pz());
            data.fillMulti<float>("pfcand_energy", packed_cand->energy());

            data.fillMulti<float>("pfcand_pt_nopuppi", packed_cand->pt());
            data.fillMulti<float>("pfcand_pt_log_nopuppi", catchInfs(std::log(packed_cand->pt()), -99));
            data.fillMulti<float>("pfcand_e_log_nopuppi", catchInfs(std::log(packed_cand->energy()), -99));

            data.fillMulti<float>("pfcand_phirel", reco::deltaPhi(*packed_cand, jet));
            data.fillMulti<float>("pfcand_etarel", etasign * (packed_cand->eta() - jet.eta()));
            data.fillMulti<float>("pfcand_abseta", std::abs(packed_cand->eta()));
            data.fillMulti<float>("pfcand_ptWrtJet", ROOT::Math::VectorUtil::Perp(packed_cand->p4().Vect(), jet.p4().Vect()));

            double x_LBGS = catchInfs(v_consti_4mom_tr.at(iConsti).py() / v_consti_4mom_tr.at(iConsti).e());
            //double x_LBGS = v_consti_4mom_tr.at(iConsti).py();
            double y_LBGS = catchInfs(v_consti_4mom_tr.at(iConsti).pz() / v_consti_4mom_tr.at(iConsti).e());
            double enFrac = catchInfs(v_consti_4mom_tr.at(iConsti).e() / jetLorentzBoost_e0);

            enFrac = (v_consti_4mom_tr.size() == 1) ? 1.0 : enFrac;
            enFrac = std::min(enFrac, 1.0);

            data.fillMulti<float>("pfcand_LBGS_x", x_LBGS);
            data.fillMulti<float>("pfcand_LBGS_y", y_LBGS);
            data.fillMulti<float>("pfcand_LBGS_enFrac", enFrac);

            data.fillMulti<float>("pfcand_puppiw", jet_helper.getPuppiWeight(cand));

            data.fillMulti<float>("pfcand_charge", packed_cand->charge());
            data.fillMulti<float>("pfcand_isEl", std::abs(packed_cand->pdgId()) == 11);
            data.fillMulti<float>("pfcand_isMu", std::abs(packed_cand->pdgId()) == 13);
            data.fillMulti<float>("pfcand_isChargedHad", std::abs(packed_cand->pdgId()) == 211);
            data.fillMulti<float>("pfcand_isGamma", std::abs(packed_cand->pdgId()) == 22);
            data.fillMulti<float>("pfcand_isNeutralHad", std::abs(packed_cand->pdgId()) == 130);

            // for neutral
            float hcal_fraction = 0.;
            if (packed_cand->pdgId() == 1 || packed_cand->pdgId() == 130)
            {
                hcal_fraction = packed_cand->hcalFraction();
            }
            else if (packed_cand->isIsolatedChargedHadron())
            {
                hcal_fraction = packed_cand->rawHcalFraction();
            }
            data.fillMulti<float>("pfcand_hcalFrac", hcal_fraction);
            data.fillMulti<float>("pfcand_hcalFracCalib", packed_cand->hcalFraction());

            // for charged
            data.fillMulti<float>("pfcand_VTX_ass", packed_cand->pvAssociationQuality());
            data.fillMulti<float>("pfcand_fromPV", packed_cand->fromPV());
            data.fillMulti<float>("pfcand_lostInnerHits", packed_cand->lostInnerHits());
            data.fillMulti<float>("pfcand_trackHighPurity", packed_cand->trackHighPurity());

            // impact parameters
            data.fillMulti<float>("pfcand_dz", catchInfs(packed_cand->dz()));
            data.fillMulti<float>("pfcand_dzsig",
                                  packed_cand->bestTrack() ? catchInfs(packed_cand->dz() / packed_cand->dzError()) : 0);
            data.fillMulti<float>("pfcand_dxy", catchInfs(packed_cand->dxy()));
            data.fillMulti<float>("pfcand_dxysig",
                                  packed_cand->bestTrack() ? catchInfs(packed_cand->dxy() / packed_cand->dxyError()) : 0);

            if (packed_cand->bestTrack())
            {
                const auto *trk = packed_cand->bestTrack();
                data.fillMulti<float>("pfcand_normchi2", catchInfs(trk->normalizedChi2()));
                data.fillMulti<float>("pfcand_quality", trk->qualityMask());
            }
            else
            {
                data.fillMulti<float>("pfcand_normchi2", 999);
                data.fillMulti<float>("pfcand_quality", 0);
            }

            // build track info map
            TrackInfoBuilder trkinfo;
            trkinfo.buildTrackInfo(builder_, *packed_cand, jet, vertices->at(0));

            // data.fillMulti<float>("pfcand_btagMomentum", catchInfs(trkinfo.getTrackMomentum()));
            // data.fillMulti<float>("pfcand_btagEta", catchInfs(trkinfo.getTrackEta()));
            data.fillMulti<float>("pfcand_btagEtaRel", catchInfs(trkinfo.getTrackEtaRel()));
            data.fillMulti<float>("pfcand_btagPtRel", catchInfs(trkinfo.getTrackPtRel()));
            // data.fillMulti<float>("pfcand_btagPPar", catchInfs(trkinfo.getTrackPPar()));
            // data.fillMulti<float>("pfcand_btagDeltaR", catchInfs(trkinfo.getTrackDeltaR()));
            data.fillMulti<float>("pfcand_btagPtRatio", catchInfs(trkinfo.getTrackPtRatio()));
            data.fillMulti<float>("pfcand_btagPParRatio", catchInfs(trkinfo.getTrackPParRatio()));
            data.fillMulti<float>("pfcand_btagSip3dVal", catchInfs(trkinfo.getTrackSip3dVal()));
            data.fillMulti<float>("pfcand_btagSip3dSig", catchInfs(trkinfo.getTrackSip3dSig()));
            data.fillMulti<float>("pfcand_btagJetDistVal", catchInfs(trkinfo.getTrackJetDistVal()));
            //    data.fillMulti<float>("pfcand_btagJetDistSig", catchInfs(trkinfo.getTrackJetDistSig()));

            // Get matching electron
            int iEl = -1;
            int nearestElIdx = -1;
            //double nearestElPt = -1;
            double nearestElDpt = 9999;
            double nearestElDR = 9999;

            for (const auto &el : electrons)
            {
                iEl++;
                double dR = ROOT::Math::VectorUtil::DeltaR(el.p4(), packed_cand->p4());
                double dpt = std::fabs(el.pt() - packed_cand->pt());
                
                //if (dR < 0.01 && el.pt() > nearestElPt)
                if (dR < 0.01 && dpt < nearestElDpt)
                {
                    nearestElIdx = iEl;
                    //nearestElPt = el.pt();
                    nearestElDpt = dpt;
                    nearestElDR = dR;
                }
            }
            
            // Get matching muon
            int iMu = -1;
            int nearestMuIdx = -1;
            //double nearestMuPt = -1;
            double nearestMuDpt = 9999;
            double nearestMuDR = 9999;

            for (const auto &mu : muons)
            {
                iMu++;
                
                if (!((mu.isGlobalMuon() || mu.isTrackerMuon()) && mu.isPFMuon()))
                {
                    continue;
                }
                
                double dR = ROOT::Math::VectorUtil::DeltaR(mu.p4(), packed_cand->p4());
                double dpt = std::fabs(mu.pt() - packed_cand->pt());

                //if (dR < 0.01 && mu.pt() > nearestMuPt)
                if (dR < 0.01 && dpt < nearestMuDpt)
                {
                    nearestMuIdx = iMu;
                    //nearestMuPt = mu.pt();
                    nearestMuDpt = dpt;
                    nearestMuDR = dR;
                }
            }
            
            // Disambiguation
            if (nearestElIdx >= 0 && nearestMuIdx >= 0)
            {
                if (std::abs(packed_cand->pdgId()) == 11)
                {
                    nearestMuIdx = -1;
                }
                
                else if (std::abs(packed_cand->pdgId()) == 13)
                {
                    nearestElIdx = -1;
                }
                
                else
                {
                    if (nearestElDR < nearestMuDR)
                    {
                        nearestMuIdx = -1;
                    }
                    
                    else
                    {
                        nearestElIdx = -1;
                    }
                }
            }
            
            data.fillMulti<float>("pfcand_isElMatched", (int) (nearestElIdx >= 0));
            data.fillMulti<float>("pfcand_isMuMatched", (int) (nearestMuIdx >= 0));
            
            // Fill electron branches
            if (nearestElIdx >= 0)
            {
                edm::Ptr<pat::Electron> elPtr(electronsH, nearestElIdx);

                for(int iVar = 0; iVar < elMvaVarManager_->getNVars(); iVar++)
                {
                    std::string varName = elMvaVarManager_->getName(iVar);
                    double val = catchInfs(elMvaVarManager_->getValue(iVar, *elPtr, extraVariables_));
                    data.fillMulti<float>(m_elVarBranchName_.at(varName), val);
                }
            }
            
            else
            {
                for(int iVar = 0; iVar < elMvaVarManager_->getNVars(); iVar++)
                {
                    std::string varName = elMvaVarManager_->getName(iVar);
                    data.fillMulti<float>(m_elVarBranchName_.at(varName), -9999.0);
                }
            }
            
            // Fill muon branches
            if (nearestMuIdx >= 0)
            {
                edm::Ptr<pat::Muon> muPtr(muonsH, nearestMuIdx);

                for(int iVar = 0; iVar < muMvaVarManager_->getNVars(); iVar++)
                {
                    std::string varName = muMvaVarManager_->getName(iVar);
                    double val = catchInfs(muMvaVarManager_->getValue(iVar, *muPtr, extraVariables_));
                    data.fillMulti<float>(m_muVarBranchName_.at(varName), val);
                }
            }
            
            else
            {
                for(int iVar = 0; iVar < muMvaVarManager_->getNVars(); iVar++)
                {
                    std::string varName = muMvaVarManager_->getName(iVar);
                    data.fillMulti<float>(m_muVarBranchName_.at(varName), -9999.0);
                }
            }
        }

        return true;
    }

} /* namespace deepntuples */
