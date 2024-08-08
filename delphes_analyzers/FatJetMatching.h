#ifndef FatJetMatching_h
#define FatJetMatching_h

#include <iostream>
#include <cassert>
#include <unordered_set>
#include <utility>
#include "TClonesArray.h"
#include "Math/LorentzVector.h"
#include <Math/Vector4D.h>
#include "classes/DelphesClasses.h"

#include "ParticleID.h"
#include "ParticleInfo.h"

class FatJetMatching {

public:
    struct FatJetMatchingResult {
        std::string label;
        std::vector<const GenParticle*> resParticles;
        std::vector<const GenParticle*> decayParticles;
        std::vector<const GenParticle*> tauDecayParticles;
        std::vector<const GenParticle*> qcdPartons;
    };
    FatJetMatchingResult& getResult() { return result_; }
    void clearResult() {
        result_.label = "Invalid";
        result_.resParticles.clear();
        result_.decayParticles.clear();
        result_.tauDecayParticles.clear();
        result_.qcdPartons.clear();
    }

public:
    FatJetMatching() {}
    FatJetMatching(double jetR, bool assignQCDLabel, bool debug) : jetR_(jetR), assignQCDLabel_(assignQCDLabel), debug_(debug) {}

    virtual ~FatJetMatching() {}

    void getLabel(const Jet *jet, const TClonesArray *branchParticle) {

        genParticles_.clear();
        for (Int_t i = 0; i < branchParticle->GetEntriesFast(); ++i) {
        genParticles_.push_back((GenParticle *)branchParticle->At(i));
        }
        processed_.clear();

        if (debug_) {
        std::cout << "\n=======\nJet (energy, pT, eta, phi) = " << jet->P4().Energy() << ", " << jet->PT << ", "
                    << jet->Eta << ", " << jet->Phi << std::endl
                    << std::endl;
        //   printGenInfoHeader();
        //   for (unsigned ipart = 0; ipart < genParticles_.size(); ++ipart) {
        //     printGenParticleInfo(genParticles_[ipart], ipart);
        //   }
        }

        for (const auto *gp : genParticles_) {
        if (processed_.count(gp))
            continue;
        processed_.insert(gp);

        auto pdgid = std::abs(gp->PID);
        if (pdgid == ParticleID::p_h0 || pdgid == ParticleID::p_H0 || pdgid == ParticleID::p_Hplus) {
            clearResult();
            res2PLabel(jet, gp);
            if (getResult().label != "Invalid") {
                return;
            }
        }
        }

        if (genParticles_.size() != processed_.size())
            throw std::logic_error("[FatJetMatching] Not all genParticles are processed!");

        if (assignQCDLabel_) {
            clearResult();
            qcdLabel(jet);
        }
    }

    int findLabelIndex() {

        if (getResult().label == "Invalid") {
            throw std::logic_error("[FatJetInfoFiller::fill]: label is Invalid");
        }
        
        int label_index = -1;
        auto it = std::find(labels_.begin(), labels_.end(), getResult().label);
        if (it != labels_.end()) {
            label_index = std::distance(labels_.begin(), it);
        }else {
            throw std::logic_error("[FatJetInfoFiller::fill]: unexpected label " + getResult().label);
        }
        return label_index;
    }

private:

    void res2PLabel(const Jet *jet, const GenParticle *parton) {

        auto res = getFinal(parton);
        getResult().resParticles.push_back(res);

        if (debug_){
            using namespace std;
            cout << "jet: (" << jet->PT << ", " << jet->Eta << ", " << jet->Phi << ", " << jet->P4().Energy() << ")" << endl;
            std::cout << "H:     "; printGenParticleInfo(res, -1);
        }

        enum XDecay {X_2p, X_tautau, X_YY, X_null};
        XDecay xdecay = X_null;
        auto daus = getDaughters(res);
        if (res->D2 - res->D1 + 1 >= 3) {
            throw std::runtime_error("[FatJetMatching::res2PLabel] X decays to >2 objects: not valid");
        }else {
            auto pdgid1 = std::abs(daus.at(0)->PID), pdgid2 = std::abs(daus.at(1)->PID);
            if ((pdgid1 == ParticleID::p_Hplus && pdgid2 == ParticleID::p_Hplus) || (pdgid1 == ParticleID::p_H0 && pdgid2 == ParticleID::p_H0)) {
                xdecay = X_YY;
            }else if (pdgid1 == ParticleID::p_tauminus && pdgid2 == ParticleID::p_tauminus) {
                xdecay = X_tautau;
            }else {
                xdecay = X_2p;
            }
        }

        if (xdecay == X_YY){
            // h->H0H0 or H+H-
            enum YMode {Y_had, Y_lep, Y_null};
            YMode yy_modes[2] = {Y_null, Y_null};
            std::vector<const GenParticle*> yy_daus;
            // found daughters of YY, and determine the their decay (had or lep),
            // then switch order to *make sure had daughters go first*
            for (int idau = res->D1; idau <= res->D2; ++idau) {
                int i = idau - res->D1;
                const auto *dau = genParticles_.at(idau);
                auto daufinal = getFinal(dau);
                getResult().resParticles.push_back(daufinal); // push the first daughter to the list (H0 or H+-)

                for (int jdau = daufinal->D1; jdau <= daufinal->D2; ++jdau){
                    const auto *ddau = genParticles_.at(jdau);
                    // determine the Y decay mode
                    if (jdau == daufinal->D1) {
                        auto dpdgid = std::abs(ddau->PID);
                        if ((dpdgid >= ParticleID::p_d && dpdgid <= ParticleID::p_b) || dpdgid == ParticleID::p_g) {
                            yy_modes[i] = Y_had;
                        }else {
                            yy_modes[i] = Y_lep;
                        }
                    }
                    yy_daus.push_back(ddau);
                }
            }
            if (yy_modes[0] == Y_lep && yy_modes[1] == Y_lep) {
                // require not both Y are leptonic
                if (debug_)
                    std::cout << "Both Ys decay leptonically, skip" << std::endl;
                return;
            }
            else if (yy_modes[0] == Y_lep && yy_modes[1] == Y_had) {
                // hadronic Y always goes first
                std::swap(yy_daus.at(0), yy_daus.at(2));
                std::swap(yy_daus.at(1), yy_daus.at(3));
                yy_modes[0] = Y_had;
                yy_modes[1] = Y_lep;
            }

            // let e/mu/tau appears before neutrinos
            if (yy_modes[1] == Y_lep) {
                auto pdgid = std::abs(yy_daus.at(2)->PID);
                if (pdgid == ParticleID::p_nu_e || pdgid == ParticleID::p_nu_mu || pdgid == ParticleID::p_nu_tau) {
                    std::swap(yy_daus.at(2), yy_daus.at(3));
                }
            }

            // go to dedicated X_YY
            // result.resParticles is completed; result.label and result.decayParticles to be assigned inside those functions
            res34PLabel(jet, yy_daus);
            return;

        }else if (xdecay == X_2p) {
            // the 2-prong case

            if (daus.size() < 2)
                throw std::logic_error("[FatJetMatching::res2PLabel] Higgs decay has less than 2 quarks!");

            double dr_q1 = deltaR(jet, daus.at(0));
            double dr_q2 = deltaR(jet, daus.at(1));
            if (dr_q1 > dr_q2){
                // swap q1 and q2 so that dr_q1<=dr_q2
                std::swap(dr_q1, dr_q2);
                std::swap(daus.at(0), daus.at(1));
            }
            auto pdgid_q1 = std::abs(daus.at(0)->PID);
            auto pdgid_q2 = std::abs(daus.at(1)->PID);

            if (debug_){
                using namespace std;
                cout << "deltaR(jet, q1)        : " << dr_q1 << endl;
                cout << "deltaR(jet, q2)        : " << dr_q2 << endl;
                cout << "pdgid(q1)              : " << pdgid_q1 << endl;
                cout << "pdgid(q2)              : " << pdgid_q2 << endl;
            }

            if (dr_q1 < jetR_ && dr_q2 < jetR_){
                getResult().decayParticles.push_back(daus.at(0));
                getResult().decayParticles.push_back(daus.at(1));

                if (pdgid_q1 == ParticleID::p_b && pdgid_q2 == ParticleID::p_b) {
                    getResult().label = "X_bb";
                }
                else if (pdgid_q1 == ParticleID::p_c && pdgid_q2 == ParticleID::p_c) {
                    getResult().label = "X_cc";
                }
                else if (pdgid_q1 == ParticleID::p_s && pdgid_q2 == ParticleID::p_s) {
                    getResult().label = "X_ss";
                }
                else if ((pdgid_q1 == ParticleID::p_u || pdgid_q1 == ParticleID::p_d) && (pdgid_q2 == ParticleID::p_u || pdgid_q2 == ParticleID::p_d)) {
                    getResult().label = "X_qq";
                }
                else if ((pdgid_q1 == ParticleID::p_b && pdgid_q2 == ParticleID::p_c) || (pdgid_q1 == ParticleID::p_c && pdgid_q2 == ParticleID::p_b)) {
                    getResult().label = "X_bc";
                }
                else if ((pdgid_q1 == ParticleID::p_b && pdgid_q2 == ParticleID::p_s) || (pdgid_q1 == ParticleID::p_s && pdgid_q2 == ParticleID::p_b)) {
                    getResult().label = "X_bs"; // this should not happen in the main sample
                }
                else if ((pdgid_q1 == ParticleID::p_c && pdgid_q2 == ParticleID::p_s) || (pdgid_q1 == ParticleID::p_s && pdgid_q2 == ParticleID::p_c)) {
                    getResult().label = "X_cs";
                }
                else if ((pdgid_q1 == ParticleID::p_b && (pdgid_q2 == ParticleID::p_u || pdgid_q2 == ParticleID::p_d)) || ((pdgid_q1 == ParticleID::p_u || pdgid_q1 == ParticleID::p_d) && pdgid_q2 == ParticleID::p_b)) {
                    getResult().label = "X_bq";
                }
                else if ((pdgid_q1 == ParticleID::p_c && (pdgid_q2 == ParticleID::p_u || pdgid_q2 == ParticleID::p_d)) || ((pdgid_q1 == ParticleID::p_u || pdgid_q1 == ParticleID::p_d) && pdgid_q2 == ParticleID::p_c)) {
                    getResult().label = "X_cq";
                }
                else if ((pdgid_q1 == ParticleID::p_s && (pdgid_q2 == ParticleID::p_u || pdgid_q2 == ParticleID::p_d)) || ((pdgid_q1 == ParticleID::p_u || pdgid_q1 == ParticleID::p_d) && pdgid_q2 == ParticleID::p_s)) {
                    getResult().label = "X_sq";
                }
                else if (pdgid_q1 == ParticleID::p_g && pdgid_q2 == ParticleID::p_g) {
                    getResult().label = "X_gg";
                }else if (pdgid_q1 == ParticleID::p_eminus && pdgid_q2 == ParticleID::p_eminus) {
                    getResult().label = "X_ee";
                }else if (pdgid_q1 == ParticleID::p_muminus && pdgid_q2 == ParticleID::p_muminus) {
                    getResult().label = "X_mm";
                }
            }
            return;

        }else if (xdecay == X_tautau) {

            if (daus.size() < 2)
                throw std::logic_error("[FatJetMatching::res2PLabel] Higgs decay has less than 2 taus!");

            // res -> tautau
            double dr_tau1 = deltaR(jet, daus.at(0));
            double dr_tau2 = deltaR(jet, daus.at(1));

            // daus_info: (daughter list, tau decay mode (0: ev, 1: mv, 2: had))
            auto tau1_daus_info = getTauDaughters(daus.at(0));
            auto tau2_daus_info = getTauDaughters(daus.at(1));

            if (debug_){
                using namespace std;
                cout << "tau1 decay ID: " << tau1_daus_info.second << endl;
                cout << "tau2 decay ID: " << tau2_daus_info.second << endl;
                cout << "deltaR(jet, tau1)        : " << dr_tau1 << endl;
                cout << "deltaR(jet, tau1-dau): " << deltaR(tau1_daus_info.first.at(0), jet) << endl;
                cout << "deltaR(jet, tau2)        : " << dr_tau2 << endl;
                cout << "deltaR(jet, tau2-dau): " << deltaR(tau2_daus_info.first.at(0), jet) << endl;
            }

            // let hadronic tau be the first one
            if (tau1_daus_info.second < 2 && tau2_daus_info.second == 2){
                std::swap(dr_tau1, dr_tau2);
                std::swap(daus.at(0), daus.at(1));
                std::swap(tau1_daus_info, tau2_daus_info);
            }
            // the first tau must be hadronic
            if (tau1_daus_info.second == 2 && dr_tau1 < jetR_){
                // push the first tau_h
                getResult().decayParticles.push_back(daus.at(0));
                // push the tau_h decay products
                getResult().tauDecayParticles.insert(getResult().tauDecayParticles.end(), tau1_daus_info.first.begin(), tau1_daus_info.first.end());

                // inspect the second tau
                if ((tau2_daus_info.second == 0 || tau2_daus_info.second == 1) && deltaR(tau2_daus_info.first.at(0), jet) < jetR_){
                    getResult().decayParticles.push_back(daus.at(1));
                    getResult().tauDecayParticles.insert(getResult().tauDecayParticles.end(), tau2_daus_info.first.begin(), tau2_daus_info.first.end());
                    if (tau2_daus_info.second == 0) {
                        getResult().label = "X_tauhtaue";
                    }else {
                        getResult().label = "X_tauhtaum";
                    }
                }else if (tau2_daus_info.second == 2 && dr_tau2 < jetR_){
                    getResult().decayParticles.push_back(daus.at(1));
                    getResult().tauDecayParticles.insert(getResult().tauDecayParticles.end(), tau2_daus_info.first.begin(), tau2_daus_info.first.end());
                    getResult().label = "X_tauhtauh";
                }
            }
            return;

        }
    }

    void res34PLabel(const Jet* jet, std::vector<const GenParticle*>& yy_daughters)
    {
        enum XDecay {X_bb, X_cc, X_ss, X_qq, X_bc, /*X_bs*,*/ X_cs, X_bq, X_cq, X_sq, X_gg, X_ee, X_mm, X_ev, X_mv, X_tauhtaue, X_tauhtaum, X_tauhtauh, X_tauev, X_taumv, X_tauhv, X_null};

        // determine the decay mode of two Ys (can be H0/H+- in defination)
        XDecay yydecay[2] = {X_null, X_null};
        int daus_matched[4] = {0, 0, 0, 0};
        getResult().decayParticles.insert(getResult().decayParticles.end(), yy_daughters.begin(), yy_daughters.end());

        for (int i = 0; i < 2; ++i) {
            const GenParticle* daus[2] = {yy_daughters.at(i*2), yy_daughters.at(i*2+1)};

            int pdgid_q1 = std::abs(daus[0]->PID);
            int pdgid_q2 = std::abs(daus[1]->PID);

            // non-tau cases
            if ((pdgid_q1 >= ParticleID::p_d && pdgid_q1 <= ParticleID::p_nu_mu) || pdgid_q1 == ParticleID::p_g) {
                
                if (pdgid_q1 == ParticleID::p_b && pdgid_q2 == ParticleID::p_b) {
                    yydecay[i] = X_bb;
                }
                else if (pdgid_q1 == ParticleID::p_c && pdgid_q2 == ParticleID::p_c) {
                    yydecay[i] = X_cc;
                }
                else if (pdgid_q1 == ParticleID::p_s && pdgid_q2 == ParticleID::p_s) {
                    yydecay[i] = X_ss;
                }
                else if ((pdgid_q1 == ParticleID::p_u || pdgid_q1 == ParticleID::p_d) && (pdgid_q2 == ParticleID::p_u || pdgid_q2 == ParticleID::p_d)) {
                    yydecay[i] = X_qq;
                }
                else if ((pdgid_q1 == ParticleID::p_b && pdgid_q2 == ParticleID::p_c) || (pdgid_q1 == ParticleID::p_c && pdgid_q2 == ParticleID::p_b)) {
                    yydecay[i] = X_bc;
                }
                // else if ((pdgid_q1 == ParticleID::p_b && pdgid_q2 == ParticleID::p_s) || (pdgid_q1 == ParticleID::p_s && pdgid_q2 == ParticleID::p_b)) {
                //     yydecay[i] = X_bs; // this should not happen
                // }
                else if ((pdgid_q1 == ParticleID::p_c && pdgid_q2 == ParticleID::p_s) || (pdgid_q1 == ParticleID::p_s && pdgid_q2 == ParticleID::p_c)) {
                    yydecay[i] = X_cs;
                }
                else if ((pdgid_q1 == ParticleID::p_b && (pdgid_q2 == ParticleID::p_u || pdgid_q2 == ParticleID::p_d)) || ((pdgid_q1 == ParticleID::p_u || pdgid_q1 == ParticleID::p_d) && pdgid_q2 == ParticleID::p_b)) {
                    yydecay[i] = X_bq;
                }
                else if ((pdgid_q1 == ParticleID::p_c && (pdgid_q2 == ParticleID::p_u || pdgid_q2 == ParticleID::p_d)) || ((pdgid_q1 == ParticleID::p_u || pdgid_q1 == ParticleID::p_d) && pdgid_q2 == ParticleID::p_c)) {
                    yydecay[i] = X_cq;
                }
                else if ((pdgid_q1 == ParticleID::p_s && (pdgid_q2 == ParticleID::p_u || pdgid_q2 == ParticleID::p_d)) || ((pdgid_q1 == ParticleID::p_u || pdgid_q1 == ParticleID::p_d) && pdgid_q2 == ParticleID::p_s)) {
                    yydecay[i] = X_sq;
                }
                else if (pdgid_q1 == ParticleID::p_g && pdgid_q2 == ParticleID::p_g) {
                    yydecay[i] = X_gg;
                }else if (pdgid_q1 == ParticleID::p_eminus && pdgid_q2 == ParticleID::p_eminus) {
                    yydecay[i] = X_ee;
                }else if (pdgid_q1 == ParticleID::p_muminus && pdgid_q2 == ParticleID::p_muminus) {
                    yydecay[i] = X_mm;
                }else if (pdgid_q1 == ParticleID::p_eminus && pdgid_q2 == ParticleID::p_nu_e) {
                    yydecay[i] = X_ev;
                }else if (pdgid_q1 == ParticleID::p_muminus && pdgid_q2 == ParticleID::p_nu_mu) {
                    yydecay[i] = X_mv;
                }

                // check if the Y daughters are matched to the jet
                if (deltaR(jet, daus[0]) < jetR_) {
                    daus_matched[i*2] = 1;
                }
                if ((pdgid_q2 != ParticleID::p_nu_e && pdgid_q2 != ParticleID::p_nu_mu) && deltaR(jet, daus[1]) < jetR_) { // ignore neutrinos: they should have matched=0
                    daus_matched[i*2+1] = 1;
                }

            }else if (pdgid_q1 == ParticleID::p_tauminus && pdgid_q2 == ParticleID::p_tauminus) {

                // H->tautau
                auto tau1_daus_info = getTauDaughters(daus[0]);
                auto tau2_daus_info = getTauDaughters(daus[1]);
                getResult().tauDecayParticles.insert(getResult().tauDecayParticles.end(), tau1_daus_info.first.begin(), tau1_daus_info.first.end());
                getResult().tauDecayParticles.insert(getResult().tauDecayParticles.end(), tau2_daus_info.first.begin(), tau2_daus_info.first.end());

                // let hadronic tau be the first one
                if (tau1_daus_info.second < 2 && tau2_daus_info.second == 2){
                    std::swap(yy_daughters.at(i*2), yy_daughters.at(i*2+1));
                    std::swap(daus[0], daus[1]);
                    std::swap(tau1_daus_info, tau2_daus_info);
                }
                if (tau1_daus_info.second == 2){ // first tau must be tau_h
                    if (tau2_daus_info.second == 0 || tau2_daus_info.second == 1){
                        if (tau2_daus_info.second == 0) {
                            yydecay[i] = X_tauhtaue;
                        }else {
                            yydecay[i] = X_tauhtaum;
                        }
                        // matching of the first tau_h and the second tau_e/tau_m
                        if (deltaR(jet, daus[0]) < jetR_) {
                            daus_matched[i*2] = 1;
                        }
                        if (deltaR(jet, tau2_daus_info.first.at(0)) < jetR_) {
                            daus_matched[i*2+1] = 1;
                        }
                    }else if (tau2_daus_info.second == 2){
                        yydecay[i] = X_tauhtauh;
                        // matching of both tau_h
                        if (deltaR(jet, daus[0]) < jetR_) {
                            daus_matched[i*2] = 1;
                        }
                        if (deltaR(jet, daus[1]) < jetR_) {
                            daus_matched[i*2+1] = 1;
                        }
                    }
                }

            }else if (pdgid_q1 == ParticleID::p_tauminus && pdgid_q2 == ParticleID::p_nu_tau) {

                // H->tauv
                auto tau_daus_info = getTauDaughters(daus[0]);
                getResult().tauDecayParticles.insert(getResult().tauDecayParticles.end(), tau_daus_info.first.begin(), tau_daus_info.first.end());

                if (tau_daus_info.second == 0 || tau_daus_info.second == 1){
                    if (tau_daus_info.second == 0) {
                        yydecay[i] = X_tauev;
                    }else {
                        yydecay[i] = X_taumv;
                    }
                    // matching of the tau_e and tau_m, neutrinos considered unmatched
                    if (deltaR(jet, tau_daus_info.first.at(0)) < jetR_) {
                        daus_matched[i*2] = 1;
                    }
                } else if (tau_daus_info.second == 2){
                    yydecay[i] = X_tauhv;
                    // matching of the tau_h
                    if (deltaR(jet, daus[0]) < jetR_) {
                        daus_matched[i*2] = 1;
                    }
                }
            }

            if (debug_){
                using namespace std;
                cout << "Y_" << i+1 << " decay mode    :" << yydecay[i] << endl;
                cout <<    "    pdgid(dau1)     : " << daus[0]->PID << endl;
                cout <<    "    pdgid(dau2)     : " << daus[1]->PID << endl;
                cout <<    "    deltaR(jet, dau1)    : " << deltaR(jet, daus[0]) << "    dau/tau-dau is matched: " << daus_matched[i*2] << endl;
                cout <<    "    deltaR(jet, dau2)    : " << deltaR(jet, daus[1]) << "    dau/tau-dau is matched: " << daus_matched[i*2+1] << endl;
                int jmax = 0;
                if (yydecay[i] == X_tauhtaue || yydecay[i] == X_tauhtaum || yydecay[i] == X_tauhtauh)
                    jmax = 2;
                else if (yydecay[i] == X_tauev || yydecay[i] == X_taumv || yydecay[i] == X_tauhv)
                    jmax = 1;
                for (int j = 0; j < jmax; j++) {
                    auto tau_daus_info = getTauDaughters(daus[j]);
                    auto tau_daus = tau_daus_info.first;
                    int tau_decay = tau_daus_info.second;
                    cout << "    tau_" << j+1 << " :" << endl;
                    cout << "        tau decay mode: " << tau_decay << endl;
                    cout << "        tau daughters: " << endl;
                    for (auto dau : tau_daus) {
                        cout << "            pdgid: " << dau->PID << "    deltaR(jet, dau): " << deltaR(jet, dau) << endl;
                    }
                }
            }

        }

        // selection criteria
        if (yydecay[0] == X_null || yydecay[1] == X_null) {
            // other Y modes will not be handled (e.g Y->tauetaue)
            return;
        }
        if (yydecay[1] == X_tauhtaue || yydecay[1] == X_tauhtaum || yydecay[1] == X_tauhtauh) {
            // requires both taus should be matched
            if (!daus_matched[2] || !daus_matched[3]) {
                return;
            }
        }
        if (daus_matched[0] + daus_matched[1] + daus_matched[2] + daus_matched[3] < 3) {
            // should have at least has three matched daughters to desinate the X_YY label
            return;
        }

        // make labels
        std::string matched_parts_str = "";
        for (int i = 0; i < 2; ++i) {
            if (yydecay[i] == X_bb || yydecay[i] == X_cc || yydecay[i] == X_ss || yydecay[i] == X_qq || 
                yydecay[i] == X_bc || yydecay[i] == X_cs || yydecay[i] == X_bq || yydecay[i] == X_cq || yydecay[i] == X_sq || 
                yydecay[i] == X_gg || yydecay[i] == X_ee || yydecay[i] == X_mm || yydecay[i] == X_ev || yydecay[i] == X_mv) {
                for (int j = 0; j < 2; ++j) {
                    // only write to string if the daughter is matched (unless it is a neutrino)
                    int pid = std::abs(yy_daughters.at(i*2+j)->PID);
                    if (daus_matched[i*2 + j] || pid == ParticleID::p_nu_e || pid == ParticleID::p_nu_mu) {
                        if (pid == ParticleID::p_b)    matched_parts_str += "b";
                        else if (pid == ParticleID::p_c)    matched_parts_str += "c";
                        else if (pid == ParticleID::p_s)    matched_parts_str += "s";
                        else if (pid == ParticleID::p_u || pid == ParticleID::p_d)    matched_parts_str += "q";
                        else if (pid == ParticleID::p_g)    matched_parts_str += "g";
                        else if (pid == ParticleID::p_eminus)    matched_parts_str += "e";
                        else if (pid == ParticleID::p_muminus)    matched_parts_str += "m";
                        else if (pid == ParticleID::p_nu_e || pid == ParticleID::p_nu_mu)    matched_parts_str += "v";
                    }
                }
            }
            else if (yydecay[i] == X_tauhtaue || yydecay[i] == X_tauhtaum || yydecay[i] == X_tauhtauh) {
                assert(daus_matched[i*2] && daus_matched[i*2+1]);
                if (yydecay[i] == X_tauhtaue)  matched_parts_str += "tauhtaue";
                else if (yydecay[i] == X_tauhtaum)    matched_parts_str += "tauhtaum";
                else if (yydecay[i] == X_tauhtauh)    matched_parts_str += "tauhtauh";
            }
            else if (yydecay[i] == X_tauev || yydecay[i] == X_taumv || yydecay[i] == X_tauhv) {
                assert(daus_matched[i*2]); // because the neutrino is not matched, all the other 3 daughters should be matched
                if (yydecay[i] == X_tauev)  matched_parts_str += "tauev";
                else if (yydecay[i] == X_taumv)    matched_parts_str += "taumv";
                else if (yydecay[i] == X_tauhv)    matched_parts_str += "tauhv";
            }
        }
        // desinate the label
        std::map<std::vector<std::string>, std::string> acceptable_strs_map = {

            // H0H0 cases (some can also be from H+H-)
            //   the first H0 -> QQ/ee/mm/tautau is all matched
            {{"bbbb"}, "X_YY_bbbb"}, {{"bbcc", "bcbc", "bccb", "cbbc", "cbcb", "ccbb"}, "X_YY_bbcc"}, {{"bbss", "bsbs", "bssb", "sbbs", "sbsb", "ssbb"}, "X_YY_bbss"}, {{"bbqq", "bqbq", "bqqb", "qbbq", "qbqb", "qqbb"}, "X_YY_bbqq"}, {{"bbgg", "ggbb"}, "X_YY_bbgg"}, {{"bbee", "eebb"}, "X_YY_bbee"}, {{"bbmm", "mmbb"}, "X_YY_bbmm"},
            {{"bbtauhtaue", "tauhtauebb"}, "X_YY_bbtauhtaue"}, {{"bbtauhtaum", "tauhtaumbb"}, "X_YY_bbtauhtaum"}, {{"bbtauhtauh", "tauhtauhbb"}, "X_YY_bbtauhtauh"},

            {{"bbb"}, "X_YY_bbb"}, {{"bbc", "bcb", "cbb"}, "X_YY_bbc"}, {{"bbs", "bsb", "sbb"}, "X_YY_bbs"}, {{"bbq", "bqb", "qbb"}, "X_YY_bbq"}, {{"bbg", "gbb"}, "X_YY_bbg"}, {{"bbe", "ebb"}, "X_YY_bbe"}, {{"bbm", "mbb"}, "X_YY_bbm"},

            {{"cccc"}, "X_YY_cccc"}, {{"ccss", "cscs", "cssc", "sccs", "scsc", "sscc"}, "X_YY_ccss"}, {{"ccqq", "cqcq", "cqqc", "qccq", "qcqc", "qqcc"}, "X_YY_ccqq"}, {{"ccgg", "ggcc"}, "X_YY_ccgg"}, {{"ccee", "eecc"}, "X_YY_ccee"}, {{"ccmm", "mmcc"}, "X_YY_ccmm"},
            {{"cctauhtaue", "tauhtauecc"}, "X_YY_cctauhtaue"}, {{"cctauhtaum", "tauhtaumcc"}, "X_YY_cctauhtaum"}, {{"cctauhtauh", "tauhtauhcc"}, "X_YY_cctauhtauh"},

            {{"ccb", "cbc", "bcc"}, "X_YY_ccb"}, {{"ccc"}, "X_YY_ccc"}, {{"ccs", "csc", "scc"}, "X_YY_ccs"}, {{"ccq", "cqc", "qcc"}, "X_YY_ccq"}, {{"ccg", "gcc"}, "X_YY_ccg"}, {{"cce", "ecc"}, "X_YY_cce"}, {{"ccm", "mcc"}, "X_YY_ccm"},

            {{"ssss"}, "X_YY_ssss"}, {{"ssqq", "sqsq", "sqqs", "qssq", "qsqs", "qqss"}, "X_YY_ssqq"}, {{"ssgg", "ggss"}, "X_YY_ssgg"}, {{"ssee", "eess"}, "X_YY_ssee"}, {{"ssmm", "mmss"}, "X_YY_ssmm"},
            {{"sstauhtaue", "tauhtauess"}, "X_YY_sstauhtaue"}, {{"sstauhtaum", "tauhtaumss"}, "X_YY_sstauhtaum"}, {{"sstauhtauh", "tauhtauhss"}, "X_YY_sstauhtauh"},

            {{"ssb", "sbs", "bss"}, "X_YY_ssb"}, {{"ssc", "scs", "css"}, "X_YY_ssc"}, {{"sss"}, "X_YY_sss"}, {{"ssq", "sqs", "qss"}, "X_YY_ssq"}, {{"ssg", "gss"}, "X_YY_ssg"}, {{"sse", "ess"}, "X_YY_sse"}, {{"ssm", "mss"}, "X_YY_ssm"},

            {{"qqqq"}, "X_YY_qqqq"}, {{"qqgg", "ggqq"}, "X_YY_qqgg"}, {{"qqee", "eeqq"}, "X_YY_qqee"}, {{"qqmm", "mmqq"}, "X_YY_qqmm"},
            {{"qqtauhtaue", "tauhtaueqq"}, "X_YY_qqtauhtaue"}, {{"qqtauhtaum", "tauhtaumqq"}, "X_YY_qqtauhtaum"}, {{"qqtauhtauh", "tauhtauhqq"}, "X_YY_qqtauhtauh"},

            {{"qqb", "qbq", "bqq"}, "X_YY_qqb"}, {{"qqc", "qcq", "cqq"}, "X_YY_qqc"}, {{"qqs", "qsq", "sqq"}, "X_YY_qqs"}, {{"qqq"}, "X_YY_qqq"}, {{"qqg", "gqq"}, "X_YY_qqg"}, {{"qqe", "eqq"}, "X_YY_qqe"}, {{"qqm", "mqq"}, "X_YY_qqm"},

            {{"gggg"}, "X_YY_gggg"}, {{"ggee", "eegg"}, "X_YY_ggee"}, {{"ggmm", "mmgg"}, "X_YY_ggmm"},
            {{"ggtauhtaue", "tauhtauegg"}, "X_YY_ggtauhtaue"}, {{"ggtauhtaum", "tauhtaumgg"}, "X_YY_ggtauhtaum"}, {{"ggtauhtauh", "tauhtauhgg"}, "X_YY_ggtauhtauh"},

            {{"ggb", "bgg"}, "X_YY_ggb"}, {{"ggc", "cgg"}, "X_YY_ggc"}, {{"ggs", "sgg"}, "X_YY_ggs"}, {{"ggq", "qgg"}, "X_YY_ggq"}, {{"ggg"}, "X_YY_ggg"}, {{"gge", "egg"}, "X_YY_gge"}, {{"ggm", "mgg"}, "X_YY_ggm"},

            {{"eeb", "bee"}, "X_YY_bee"}, {{"eec", "cee"}, "X_YY_cee"}, {{"ees", "see"}, "X_YY_see"}, {{"eeq", "qee"}, "X_YY_qee"}, {{"eeg", "gee"}, "X_YY_gee"},
            {{"mmb", "bmm"}, "X_YY_bmm"}, {{"mmc", "cmm"}, "X_YY_cmm"}, {{"mms", "smm"}, "X_YY_smm"}, {{"mmq", "qmm"}, "X_YY_qmm"}, {{"mmg", "gmm"}, "X_YY_gmm"},

            {{"tauhtaueb", "btauhtaue"}, "X_YY_btauhtaue"}, {{"tauhtauec", "ctauhtaue"}, "X_YY_ctauhtaue"}, {{"tauhtaues", "stauhtaue"}, "X_YY_stauhtaue"}, {{"tauhtaueq", "qtauhtaue"}, "X_YY_qtauhtaue"}, {{"tauhtaueg", "gtauhtaue"}, "X_YY_gtauhtaue"},
            {{"tauhtaumb", "btauhtaum"}, "X_YY_btauhtaum"}, {{"tauhtaumc", "ctauhtaum"}, "X_YY_ctauhtaum"}, {{"tauhtaums", "stauhtaum"}, "X_YY_stauhtaum"}, {{"tauhtaumq", "qtauhtaum"}, "X_YY_qtauhtaum"}, {{"tauhtaumg", "gtauhtaum"}, "X_YY_gtauhtaum"},
            {{"tauhtauhb", "btauhtauh"}, "X_YY_btauhtauh"}, {{"tauhtauhc", "ctauhtauh"}, "X_YY_ctauhtauh"}, {{"tauhtauhs", "stauhtauh"}, "X_YY_stauhtauh"}, {{"tauhtauhq", "qtauhtauh"}, "X_YY_qtauhtauh"}, {{"tauhtauhg", "gtauhtauh"}, "X_YY_gtauhtauh"},

            // H+H- case
            {{"qqbq", "qqqb", "bqqq", "qbqq"}, "X_YY_qqqb"}, {{"qqcq", "qqqc", "cqqq", "qcqq"}, "X_YY_qqqc"}, {{"qqsq", "qqqs", "sqqq", "qsqq"}, "X_YY_qqqs"},
            {{"bbcs", "bbsc", "bcbs", "bcsb", "bsbc", "bscb", "cbbs", "cbsb", "csbb", "sbbc", "sbcb", "scbb"}, "X_YY_bbcs"}, {{"bbcq", "bbqc", "bcbq", "bcqb", "bqbc", "bqcb", "cbbq", "cbqb", "cqbb", "qbbc", "qbcb", "qcbb"}, "X_YY_bbcq"}, {{"bbsq", "bbqs", "bsbq", "bsqb", "bqbs", "bqsb", "sbbq", "sbqb", "sqbb", "qbbs", "qbsb", "qsbb"}, "X_YY_bbsq"},
            {{"ccbs", "ccsb", "cbcs", "cbsc", "cscb", "csbc", "bccs", "bcsc", "bscc", "sccb", "scbc", "sbcc"}, "X_YY_ccbs"}, {{"ccbq", "ccqb", "cbcq", "cbqc", "cqcb", "cqbc", "bccq", "bcqc", "bqcc", "qccb", "qcbc", "qbcc"}, "X_YY_ccbq"}, {{"ccsq", "ccqs", "cscq", "csqc", "cqcs", "cqsc", "sccq", "scqc", "sqcc", "qccs", "qcsc", "qscc"}, "X_YY_ccsq"},
            {{"ssbc", "sscb", "sbsc", "sbcs", "scsb", "scbs", "bssc", "bscs", "bcss", "cssb", "csbs", "cbss"}, "X_YY_ssbc"}, {{"ssbq", "ssqb", "sbsq", "sbqs", "sqsb", "sqbs", "bssq", "bsqs", "bqss", "qssb", "qsbs", "qbss"}, "X_YY_ssbq"}, {{"sscq", "ssqc", "scsq", "scqs", "sqsc", "sqcs", "cssq", "csqs", "cqss", "qssc", "qscs", "qcss"}, "X_YY_sscq"},
            {{"qqbc", "qqcb", "qbqc", "qbcq", "qcqb", "qcbq", "bqqc", "bqcq", "bcqq", "cqqb", "cqbq", "cbqq"}, "X_YY_qqbc"}, {{"qqbs", "qqsb", "qbqs", "qbsq", "qsqb", "qsbq", "bqqs", "bqsq", "bsqq", "sqqb", "sqbq", "sbqq"}, "X_YY_qqbs"}, {{"qqcs", "qqsc", "qcqs", "qcsq", "qsqc", "qscq", "cqqs", "cqsq", "csqq", "sqqc", "sqcq", "scqq"}, "X_YY_qqcs"},
            {{"bcsq", "bcqs", "bscq", "bsqc", "bqcs", "bqsc", "cbsq", "cbqs", "csbq", "csqb", "cqbs", "cqsb", "scbq", "scqb", "sbcq", "sbqc", "sqcb", "sqbc", "qcsb", "qcbs", "qscb", "qsbc", "qbcs", "qbsc"}, "X_YY_bcsq"},

            {{"bcs", "bsc", "cbs", "csb", "sbc", "scb"}, "X_YY_bcs"}, {{"bcq", "bqc", "cbq", "cqb", "qbc", "qcb"}, "X_YY_bcq"}, {{"bsq", "bqs", "sbq", "sqb", "qbs", "qsb"}, "X_YY_bsq"}, {{"csq", "cqs", "scq", "sqc", "qcs", "qsc"}, "X_YY_csq"}, 

            {{"bcev", "cbev", "evbc", "evcb"}, "X_YY_bcev"}, {{"csev", "scev", "evcs", "evsc"}, "X_YY_csev"}, {{"bqev", "qbev", "evbq", "evqb"}, "X_YY_bqev"}, {{"cqev", "qcev", "evcq", "evqc"}, "X_YY_cqev"}, {{"sqev", "qsev", "evsq", "evqs"}, "X_YY_sqev"}, {{"qqev", "evqq"}, "X_YY_qqev"},
            {{"bcmv", "cbmv", "mvbc", "mvcb"}, "X_YY_bcmv"}, {{"csmv", "scmv", "mvcs", "mvsc"}, "X_YY_csmv"}, {{"bqmv", "qbmv", "mvbq", "mvqb"}, "X_YY_bqmv"}, {{"cqmv", "qcmv", "mvcq", "mvqc"}, "X_YY_cqmv"}, {{"sqmv", "qsmv", "mvsq", "mvqs"}, "X_YY_sqmv"}, {{"qqmv", "mvqq"}, "X_YY_qqmv"},
            {{"bctauev", "cbtauev", "tauevbc", "tauevcb"}, "X_YY_bctauev"}, {{"cstauev", "sctauev", "tauevcs", "tauevsc"}, "X_YY_cstauev"}, {{"bqtauev", "qbtauev", "tauevbq", "tauevqb"}, "X_YY_bqtauev"}, {{"cqtauev", "qctauev", "tauevcq", "tauevqc"}, "X_YY_cqtauev"}, {{"sqtauev", "qstauev", "tauevsq", "tauevqs"}, "X_YY_sqtauev"}, {{"qqtauev", "tauevqq"}, "X_YY_qqtauev"},
            {{"bctaumv", "cbtaumv", "taumvbc", "taumvcb"}, "X_YY_bctaumv"}, {{"cstaumv", "sctaumv", "taumvcs", "taumvsc"}, "X_YY_cstaumv"}, {{"bqtaumv", "qbtaumv", "taumvbq", "taumvqb"}, "X_YY_bqtaumv"}, {{"cqtaumv", "qctaumv", "taumvcq", "taumvqc"}, "X_YY_cqtaumv"}, {{"sqtaumv", "qstaumv", "taumvsq", "taumvqs"}, "X_YY_sqtaumv"}, {{"qqtaumv", "taumvqq"}, "X_YY_qqtaumv"},
            {{"bctauhv", "cbtauhv", "tauhvbc", "tauhvcb"}, "X_YY_bctauhv"}, {{"cstauhv", "sctauhv", "tauhvcs", "tauhvsc"}, "X_YY_cstauhv"}, {{"bqtauhv", "qbtauhv", "tauhvbq", "tauhvqb"}, "X_YY_bqtauhv"}, {{"cqtauhv", "qctauhv", "tauhvcq", "tauhvqc"}, "X_YY_cqtauhv"}, {{"sqtauhv", "qstauhv", "tauhvsq", "tauhvqs"}, "X_YY_sqtauhv"}, {{"qqtauhv", "tauhvqq"}, "X_YY_qqtauhv"},
        };

        for (auto it = acceptable_strs_map.begin(); it != acceptable_strs_map.end(); ++it) {
            auto& acceptable_strs = it->first;
            std::string& label = it->second;

            // iterate over all string options
            for (auto it_str = acceptable_strs.begin(); it_str != acceptable_strs.end(); ++it_str) {
                if (matched_parts_str == (*it_str)) {
                    getResult().label = label; // assign the label!
                    return;
                }
            }
        }
        throw std::logic_error("[FatJetMatching::res34PLabel] Unmatched YY label: " + matched_parts_str);
    }

    void qcdLabel(const Jet* jet) {

        int n_b=0, n_c=0, n_s=0;
        for (const auto *gp : genParticles_) {
            int pdgid = std::abs(gp->PID);
            // select quarks right before hadronization (status=71) and pT > 10
            if (pdgid >= ParticleID::p_d && pdgid <= ParticleID::p_b && gp->Status == 71 && gp->PT > 10 && deltaR(jet, gp) < jetR_) {
                
                if (debug_) {
                    std::cout << "Matched quark "; printGenParticleInfo(gp, -1);
                }
                getResult().qcdPartons.push_back(gp);
                if (pdgid == ParticleID::p_b)    n_b++;
                else if (pdgid == ParticleID::p_c)    n_c++;
                else if (pdgid == ParticleID::p_s)    n_s++;
            }
        }
        if (n_b >= 2 && n_c >= 2 && n_s >= 2)  getResult().label = "QCD_bbccss";
        else if (n_b >= 2 && n_c >= 2 && n_s == 1)  getResult().label = "QCD_bbccs";
        else if (n_b >= 2 && n_c >= 2 && n_s == 0)  getResult().label = "QCD_bbcc";
        else if (n_b >= 2 && n_c == 1 && n_s >= 2)  getResult().label = "QCD_bbcss";
        else if (n_b >= 2 && n_c == 1 && n_s == 1)  getResult().label = "QCD_bbcs";
        else if (n_b >= 2 && n_c == 1 && n_s == 0)  getResult().label = "QCD_bbc";
        else if (n_b >= 2 && n_c == 0 && n_s >= 2)  getResult().label = "QCD_bbss";
        else if (n_b >= 2 && n_c == 0 && n_s == 1)  getResult().label = "QCD_bbs";
        else if (n_b >= 2 && n_c == 0 && n_s == 0)  getResult().label = "QCD_bb";
        else if (n_b == 1 && n_c >= 2 && n_s >= 2)  getResult().label = "QCD_bccss";
        else if (n_b == 1 && n_c >= 2 && n_s == 1)  getResult().label = "QCD_bccs";
        else if (n_b == 1 && n_c >= 2 && n_s == 0)  getResult().label = "QCD_bcc";
        else if (n_b == 1 && n_c == 1 && n_s >= 2)  getResult().label = "QCD_bcss";
        else if (n_b == 1 && n_c == 1 && n_s == 1)  getResult().label = "QCD_bcs";
        else if (n_b == 1 && n_c == 1 && n_s == 0)  getResult().label = "QCD_bc";
        else if (n_b == 1 && n_c == 0 && n_s >= 2)  getResult().label = "QCD_bss";
        else if (n_b == 1 && n_c == 0 && n_s == 1)  getResult().label = "QCD_bs";
        else if (n_b == 1 && n_c == 0 && n_s == 0)  getResult().label = "QCD_b";
        else if (n_b == 0 && n_c >= 2 && n_s >= 2)  getResult().label = "QCD_ccss";
        else if (n_b == 0 && n_c >= 2 && n_s == 1)  getResult().label = "QCD_ccs";
        else if (n_b == 0 && n_c >= 2 && n_s == 0)  getResult().label = "QCD_cc";
        else if (n_b == 0 && n_c == 1 && n_s >= 2)  getResult().label = "QCD_css";
        else if (n_b == 0 && n_c == 1 && n_s == 1)  getResult().label = "QCD_cs";
        else if (n_b == 0 && n_c == 1 && n_s == 0)  getResult().label = "QCD_c";
        else if (n_b == 0 && n_c == 0 && n_s >= 2)  getResult().label = "QCD_ss";
        else if (n_b == 0 && n_c == 0 && n_s == 1)  getResult().label = "QCD_s";
        else if (n_b == 0 && n_c == 0 && n_s == 0)  getResult().label = "QCD_light";

    }

private:
    void printGenInfoHeader() const {
        using namespace std;
        cout << right << setw(6) << "#"
            << " " << setw(10) << "pdgId"
            << "  "
            << "Chg"
            << "  " << setw(10) << "Mass"
            << "  " << setw(48) << " Momentum" << left << "  " << setw(10) << "Mothers"
            << " " << setw(30) << "Daughters" << endl;
    }

    void printGenParticleInfo(const GenParticle *genParticle, const int idx) const {
        using namespace std;
        cout << right << setw(3) << genParticle->Status;
        cout << right << setw(3) << idx << " " << setw(10) << genParticle->PID << "  ";
        cout << right << "  " << setw(3) << genParticle->Charge << "  "
             << TString::Format("%10.3g", genParticle->Mass < 1e-5 ? 0 : genParticle->Mass);
        cout << left << setw(50)
             << TString::Format("  (E=%6.4g pT=%6.4g eta=%7.3g phi=%7.3g)",
                                genParticle->P4().Energy(),
                                genParticle->PT,
                                genParticle->Eta,
                                genParticle->Phi);

        TString mothers;
        if (genParticle->M1 >= 0) {
            mothers += genParticle->M1;
        }
        if (genParticle->M2 >= 0) {
            mothers += ",";
            mothers += genParticle->M2;
        }
        cout << "  " << setw(10) << mothers;

        TString daughters;
        for (int iDau = genParticle->D1; iDau <= genParticle->D2; ++iDau) {
            if (daughters.Length())
                daughters += ",";
            daughters += iDau;
        }
        cout << " " << setw(30) << daughters << endl;
    }

    const GenParticle *getFinal(const GenParticle *particle) {
        // will mark intermediate particles as processed
        if (!particle)
            return nullptr;
        processed_.insert(particle);
        const GenParticle *final = particle;

        while (final->D1 >= 0) {
            const GenParticle *chain = nullptr;
            for (int idau = final->D1; idau <= final->D2; ++idau) {
                if (genParticles_.at(idau)->PID == particle->PID) {
                    chain = genParticles_.at(idau);
                    processed_.insert(chain);
                    break;
                }
            }
            if (!chain)
                break;
            final = chain;
        }
        return final;
    }

    bool isHadronic(const GenParticle *particle, bool allow_gluon = false) const {
        // particle needs to be the final version before decay
        if (!particle)
            throw std::invalid_argument("[FatJetMatching::isHadronic()] Null particle!");
        for (const auto *dau : getDaughters(particle)) {
            auto pdgid = std::abs(dau->PID);
            if (pdgid >= ParticleID::p_d && pdgid <= ParticleID::p_b)
                return true;
            if (allow_gluon && pdgid == ParticleID::p_g)
                return true;
        }
        return false;
    }

    std::vector<const GenParticle *> getDaughters(const GenParticle *particle) const {
        std::vector<const GenParticle *> daughters;
        for (int idau = particle->D1; idau <= particle->D2; ++idau) {
            daughters.push_back(genParticles_.at(idau));
        }
        return daughters;
    }

    std::vector<const GenParticle *> getDaughterQuarks(const GenParticle *particle, bool allow_gluon = false) {
        std::vector<const GenParticle *> daughters;
        for (int idau = particle->D1; idau <= particle->D2; ++idau) {
            const auto *dau = genParticles_.at(idau);
            auto pdgid = std::abs(dau->PID);
            if (pdgid >= ParticleID::p_d && pdgid <= ParticleID::p_b) {
                daughters.push_back(dau);
            }
            if (allow_gluon && pdgid == ParticleID::p_g) {
                daughters.push_back(dau);
            }
        }
        return daughters;
    }

    std::pair<std::vector<const GenParticle*>, int> getTauDaughters(const GenParticle* particle) {
        const auto tau = getFinal(particle);
        auto daughters = getDaughters(tau);
        for (const auto & dau: daughters){
            auto pdgid = std::abs(dau->PID);
            if (pdgid == ParticleID::p_eminus)  return std::make_pair(daughters, 0);
            if (pdgid == ParticleID::p_muminus)  return std::make_pair(daughters, 1);
        }
        return std::make_pair(daughters, 2); // hadronic mode
    }


private:
    double jetR_ = 0.8;
    bool assignQCDLabel_ = false;
    bool debug_ = false;
    std::vector<const GenParticle *> genParticles_;
    std::unordered_set<const GenParticle *> processed_;
    FatJetMatchingResult result_{"Invalid", std::vector<const GenParticle*>(), std::vector<const GenParticle*>(), std::vector<const GenParticle*>(), std::vector<const GenParticle*>()};

    std::vector<std::string> labels_{
        // X->2-prong
        "X_bb", "X_cc", "X_ss", "X_qq", "X_bc", "X_cs", "X_bq", "X_cq", "X_sq", "X_gg", "X_ee", "X_mm", "X_tauhtaue", "X_tauhtaum", "X_tauhtauh",

        // X->YY->3/4-prong
        "X_YY_bbbb", "X_YY_bbcc", "X_YY_bbss", "X_YY_bbqq", "X_YY_bbgg", "X_YY_bbee", "X_YY_bbmm",
        "X_YY_bbtauhtaue", "X_YY_bbtauhtaum", "X_YY_bbtauhtauh",

        "X_YY_bbb", "X_YY_bbc", "X_YY_bbs", "X_YY_bbq", "X_YY_bbg", "X_YY_bbe", "X_YY_bbm",

        "X_YY_cccc", "X_YY_ccss", "X_YY_ccqq", "X_YY_ccgg", "X_YY_ccee", "X_YY_ccmm",
        "X_YY_cctauhtaue", "X_YY_cctauhtaum", "X_YY_cctauhtauh",

        "X_YY_ccb", "X_YY_ccc", "X_YY_ccs", "X_YY_ccq", "X_YY_ccg", "X_YY_cce", "X_YY_ccm",

        "X_YY_ssss", "X_YY_ssqq", "X_YY_ssgg", "X_YY_ssee", "X_YY_ssmm",
        "X_YY_sstauhtaue", "X_YY_sstauhtaum", "X_YY_sstauhtauh",

        "X_YY_ssb", "X_YY_ssc", "X_YY_sss", "X_YY_ssq", "X_YY_ssg", "X_YY_sse", "X_YY_ssm",

        "X_YY_qqqq", "X_YY_qqgg", "X_YY_qqee", "X_YY_qqmm",
        "X_YY_qqtauhtaue", "X_YY_qqtauhtaum", "X_YY_qqtauhtauh",

        "X_YY_qqb", "X_YY_qqc", "X_YY_qqs", "X_YY_qqq", "X_YY_qqg", "X_YY_qqe", "X_YY_qqm",

        "X_YY_gggg", "X_YY_ggee", "X_YY_ggmm",
        "X_YY_ggtauhtaue", "X_YY_ggtauhtaum", "X_YY_ggtauhtauh",

        "X_YY_ggb", "X_YY_ggc", "X_YY_ggs", "X_YY_ggq", "X_YY_ggg", "X_YY_gge", "X_YY_ggm",

        "X_YY_bee", "X_YY_cee", "X_YY_see", "X_YY_qee", "X_YY_gee",
        "X_YY_bmm", "X_YY_cmm", "X_YY_smm", "X_YY_qmm", "X_YY_gmm",

        "X_YY_btauhtaue", "X_YY_ctauhtaue", "X_YY_stauhtaue", "X_YY_qtauhtaue", "X_YY_gtauhtaue",
        "X_YY_btauhtaum", "X_YY_ctauhtaum", "X_YY_stauhtaum", "X_YY_qtauhtaum", "X_YY_gtauhtaum",
        "X_YY_btauhtauh", "X_YY_ctauhtauh", "X_YY_stauhtauh", "X_YY_qtauhtauh", "X_YY_gtauhtauh",

        // (note: for H+H- decays, X_YY_bbcs, X_YY_bbsq, X_YY_ssbc, X_YY_ssbq are not available)
        "X_YY_qqqb", "X_YY_qqqc", "X_YY_qqqs",
        "X_YY_bbcq",
        "X_YY_ccbs", "X_YY_ccbq", "X_YY_ccsq",
        "X_YY_sscq",
        "X_YY_qqbc", "X_YY_qqbs", "X_YY_qqcs",
        "X_YY_bcsq",

        "X_YY_bcs", "X_YY_bcq", "X_YY_bsq", "X_YY_csq", 

        "X_YY_bcev", "X_YY_csev", "X_YY_bqev", "X_YY_cqev", "X_YY_sqev", "X_YY_qqev",
        "X_YY_bcmv", "X_YY_csmv", "X_YY_bqmv", "X_YY_cqmv", "X_YY_sqmv", "X_YY_qqmv",
        "X_YY_bctauev", "X_YY_cstauev", "X_YY_bqtauev", "X_YY_cqtauev", "X_YY_sqtauev", "X_YY_qqtauev",
        "X_YY_bctaumv", "X_YY_cstaumv", "X_YY_bqtaumv", "X_YY_cqtaumv", "X_YY_sqtaumv", "X_YY_qqtaumv",
        "X_YY_bctauhv", "X_YY_cstauhv", "X_YY_bqtauhv", "X_YY_cqtauhv", "X_YY_sqtauhv", "X_YY_qqtauhv",

        // QCD
        "QCD_bbccss", "QCD_bbccs", "QCD_bbcc", "QCD_bbcss", "QCD_bbcs", "QCD_bbc", "QCD_bbss", "QCD_bbs", "QCD_bb",
        "QCD_bccss", "QCD_bccs", "QCD_bcc", "QCD_bcss", "QCD_bcs", "QCD_bc", "QCD_bss", "QCD_bs", "QCD_b",
        "QCD_ccss", "QCD_ccs", "QCD_cc", "QCD_css", "QCD_cs", "QCD_c", "QCD_ss", "QCD_s", "QCD_light",
    };
};

#endif