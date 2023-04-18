#include "DeepNTuples/NtupleCommons/interface/CommonUtils.h"

float commonutils::get_cosThetaStar(math::XYZTLorentzVector p4_top, math::XYZTLorentzVector p4_b, math::XYZTLorentzVector p4_dau1, math::XYZTLorentzVector p4_dau2, int topType) {
    // Only for hadronic tops
    if(topType == 0 && ROOT::Math::VectorUtil::InvariantMass(p4_b, p4_dau2) < ROOT::Math::VectorUtil::InvariantMass(p4_b, p4_dau1)){
        std::swap(p4_dau1, p4_dau2);
    }
    
    auto boostVector = p4_top.BoostToCM();
    auto p4_dau1_boosted = ROOT::Math::VectorUtil::boost(p4_dau1, boostVector);
    
    //printf("top: (%0.4f, %0.4f, %0.4f, %0.4f, %0.4f)\n", p4_top.px(), p4_top.py(), p4_top.pz(), p4_top.e(), p4_top.mass());
    
    //auto p4_top_rest = ROOT::Math::VectorUtil::boost(p4_top, boostVector);
    //printf("  |-- top_rest: (%0.4f, %0.4f, %0.4f, %0.4f, %0.4f)\n", p4_top_rest.px(), p4_top_rest.py(), p4_top_rest.pz(), p4_top_rest.e(), p4_top_rest.mass());
    
    //auto p4_sum = p4_b+p4_dau1+p4_dau2;
    //printf("  |-- b+q1+q2: (%0.4f, %0.4f, %0.4f, %0.4f, %0.4f)\n", p4_sum.px(), p4_sum.py(), p4_sum.pz(), p4_sum.e(), p4_sum.mass());
    
    float cosThetaStar = p4_top.Vect().Unit().Dot(p4_dau1_boosted.Vect().Unit());
    
    return cosThetaStar;
}