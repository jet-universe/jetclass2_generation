#include <iostream>
#include <unordered_set>
#include <utility>
#include "TClonesArray.h"
#include "classes/DelphesClasses.h"
#include "ExRootAnalysis/ExRootTreeReader.h"

#include "FatJetMatching.h"
#include "EventData.h"

// #ifdef __CLING__
// R__LOAD_LIBRARY(libDelphes)
// #include "classes/DelphesClasses.h"
// #include "external/ExRootAnalysis/ExRootTreeReader.h"
// #else
// class ExRootTreeReader;
// #endif

void makeNtuples(TString inputFile, TString outputFile, TString jetBranch = "JetPUPPIAK8", TString genjetBranch = "GenJetAK8", bool assignQCDLabel = false, bool debug = false) {
    // gSystem->Load("libDelphes");

    TFile *fout = new TFile(outputFile, "RECREATE");
    TTree *tree = new TTree("tree", "tree");

    // define branches
    std::vector<std::pair<std::string, std::string>> branchList = {
        // particle (jet constituent) features
        {"part_px", "vector<float>"},
        {"part_py", "vector<float>"},
        {"part_pz", "vector<float>"},
        {"part_energy", "vector<float>"},
        {"part_deta", "vector<float>"},
        {"part_dphi", "vector<float>"},
        {"part_d0val", "vector<float>"},
        {"part_d0err", "vector<float>"},
        {"part_dzval", "vector<float>"},
        {"part_dzerr", "vector<float>"},
        {"part_charge", "vector<int>"},
        {"part_isElectron", "vector<bool>"},
        {"part_isMuon", "vector<bool>"},
        {"part_isPhoton", "vector<bool>"},
        {"part_isChargedHadron", "vector<bool>"},
        {"part_isNeutralHadron", "vector<bool>"},
        // jet features
        {"jet_pt", "float"},
        {"jet_eta", "float"},
        {"jet_phi", "float"},
        {"jet_energy", "float"},
        {"jet_sdmass", "float"},
        {"jet_nparticles", "int"},
        {"jet_tau1", "float"},
        {"jet_tau2", "float"},
        {"jet_tau3", "float"},
        {"jet_tau4", "float"},
        {"jet_label", "int"},
        // gen-particle (genjet constituent) features
        {"genpart_px", "vector<float>"},
        {"genpart_py", "vector<float>"},
        {"genpart_pz", "vector<float>"},
        {"genpart_energy", "vector<float>"},
        {"genpart_jet_deta", "vector<float>"},
        {"genpart_jet_dphi", "vector<float>"},
        {"genpart_x", "vector<float>"},
        {"genpart_y", "vector<float>"},
        {"genpart_z", "vector<float>"},
        {"genpart_t", "vector<float>"},
        {"genpart_pid", "vector<int>"},
        // genjet features
        {"genjet_pt", "float"},
        {"genjet_eta", "float"},
        {"genjet_phi", "float"},
        {"genjet_energy", "float"},
        {"genjet_sdmass", "float"},
        {"genjet_nparticles", "int"},
        // aux genpart features
        {"aux_genpart_pt", "vector<float>"},
        {"aux_genpart_eta", "vector<float>"},
        {"aux_genpart_phi", "vector<float>"},
        {"aux_genpart_mass", "vector<float>"},
        {"aux_genpart_pid", "vector<int>"},
        {"aux_genpart_isResX", "vector<bool>"},
        {"aux_genpart_isResY", "vector<bool>"},
        {"aux_genpart_isResDecayProd", "vector<bool>"},
        {"aux_genpart_isTauDecayProd", "vector<bool>"},
        {"aux_genpart_isQcdParton", "vector<bool>"}
    };
    EventData data(branchList);
    data.setOutputBranch(tree);

    // Read input
    TChain *chain = new TChain("Delphes");
    chain->Add(inputFile);
    ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
    Long64_t allEntries = treeReader->GetEntries();

    std::cerr << "** Input file:    " << inputFile << std::endl;
    std::cerr << "** Jet branch:    " << jetBranch << std::endl;
    std::cerr << "** GenJet branch: " << genjetBranch << std::endl;
    std::cerr << "** Total events:  " << allEntries << std::endl;

    // Analyze
    TClonesArray *branchVertex = treeReader->UseBranch("Vertex"); // used for pileup
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchPFCand = treeReader->UseBranch("ParticleFlowCandidate");
    TClonesArray *branchJet = treeReader->UseBranch(jetBranch);
    TClonesArray *branchGenJet = treeReader->UseBranch(genjetBranch);

    double jetR = jetBranch.Contains("AK15") ? 1.5 : 0.8;
    std::cerr << "jetR = " << jetR << std::endl;

    FatJetMatching fjmatch(jetR, assignQCDLabel, debug);

    // Loop over all events
    int num_processed = 0;
    for (Long64_t entry = 0; entry < allEntries; ++entry) {
        if (entry % 1000 == 0) {
            std::cerr << "processing " << entry << " of " << allEntries << " events." << std::endl;
        }

        // Load selected branches with data from specified event
        treeReader->ReadEntry(entry);

        // Loop over all jets in event
        std::vector<int> genjet_used_inds = {};
        for (Int_t i = 0; i < branchJet->GetEntriesFast(); ++i) {
            const Jet *jet = (Jet *)branchJet->At(i);

            if (jet->PT < 120 || std::abs(jet->Eta) > 2.5)
                continue;

            data.reset();

            // Get the GEN label
            fjmatch.getLabel(jet, branchParticle);
            auto fjlabel = fjmatch.getResult().label;

            if (fjlabel == "Invalid")
                continue;

            if (debug) {
                std::cerr << ">> debug fjlabel     : " << fjlabel << "  " << std::endl;
                std::cerr << "   resParticles      : "; for (auto& p: fjmatch.getResult().resParticles) {std::cerr << p->PID << " ";} std::cerr << std::endl;
                std::cerr << "   decayParticles    : "; for (auto& p: fjmatch.getResult().decayParticles) {std::cerr << p->PID << " ";} std::cerr << std::endl;
                std::cerr << "   tauDecayParticles : "; for (auto& p: fjmatch.getResult().tauDecayParticles) {std::cerr << p->PID << " ";} std::cerr << std::endl;
                std::cerr << "   qcdPartons        : "; for (auto& p: fjmatch.getResult().qcdPartons) {std::cerr << p->PID << " ";} std::cerr << std::endl;
            }

            // GEN label and original particles (as auxiliary vars)
            data.intVars.at("jet_label") = fjmatch.findLabelIndex();

            auto fillAuxVars = [](EventData& data, const auto& part, bool isResX, bool isResY, bool isResDecayProd, bool isTauDecayProd, bool isQcdParton) {
                data.vfloatVars.at("aux_genpart_pt")->push_back(part->PT);
                data.vfloatVars.at("aux_genpart_eta")->push_back(part->Eta);
                data.vfloatVars.at("aux_genpart_phi")->push_back(part->Phi);
                data.vfloatVars.at("aux_genpart_mass")->push_back(part->Mass);
                data.vintVars.at("aux_genpart_pid")->push_back(part->PID);
                data.vboolVars.at("aux_genpart_isResX")->push_back(isResX);
                data.vboolVars.at("aux_genpart_isResY")->push_back(isResY);
                data.vboolVars.at("aux_genpart_isResDecayProd")->push_back(isResDecayProd);
                data.vboolVars.at("aux_genpart_isTauDecayProd")->push_back(isTauDecayProd);
                data.vboolVars.at("aux_genpart_isQcdParton")->push_back(isQcdParton);
            };
            int nRes = 0;
            for (const auto &p : fjmatch.getResult().resParticles) {
                fillAuxVars(data, p, nRes == 0, nRes > 0, false, false, false);
                ++nRes;
            }
            for (const auto &p : fjmatch.getResult().decayParticles) {
                fillAuxVars(data, p, false, false, true, false, false);
            }
            for (const auto &p : fjmatch.getResult().tauDecayParticles) {
                fillAuxVars(data, p, false, false, false, true, false);
            }
            for (const auto &p : fjmatch.getResult().qcdPartons) {
                fillAuxVars(data, p, false, false, false, false, true);
            }

            // Jet features
            data.floatVars.at("jet_pt") = jet->PT;
            data.floatVars.at("jet_eta") = jet->Eta;
            data.floatVars.at("jet_phi") = jet->Phi;
            data.floatVars.at("jet_energy") = jet->P4().Energy();

            data.floatVars.at("jet_sdmass") = jet->SoftDroppedP4[0].M();
            data.floatVars.at("jet_tau1") = jet->Tau[0];
            data.floatVars.at("jet_tau2") = jet->Tau[1];
            data.floatVars.at("jet_tau3") = jet->Tau[2];
            data.floatVars.at("jet_tau4") = jet->Tau[3];

            // Loop over all jet's constituents
            std::vector<ParticleInfo> particles;
            for (Int_t j = 0; j < jet->Constituents.GetEntriesFast(); ++j) {
                const TObject *object = jet->Constituents.At(j);

                // Check if the constituent is accessible
                if (!object)
                    continue;

                if (object->IsA() == GenParticle::Class()) {
                    particles.emplace_back((GenParticle *)object);
                } else if (object->IsA() == ParticleFlowCandidate::Class()) {
                    particles.emplace_back((ParticleFlowCandidate *)object);
                }
                const auto &p = particles.back();
                if (std::abs(p.pz) > 10000 || std::abs(p.eta) > 5 || p.pt <= 0) {
                    particles.pop_back();
                }
            }

            // sort particles by pt
            std::sort(particles.begin(), particles.end(), [](const auto &a, const auto &b) { return a.pt > b.pt; });

            // Load the primary vertex
            const Vertex *pv = (branchVertex != nullptr) ? ((Vertex *)branchVertex->At(0)) : nullptr;

            data.intVars["jet_nparticles"] = particles.size();
            for (const auto &p : particles) {
                data.vfloatVars.at("part_px")->push_back(p.px);
                data.vfloatVars.at("part_py")->push_back(p.py);
                data.vfloatVars.at("part_pz")->push_back(p.pz);
                data.vfloatVars.at("part_energy")->push_back(p.energy);
                data.vfloatVars.at("part_deta")->push_back((jet->Eta > 0 ? 1 : -1) * (p.eta - jet->Eta));
                data.vfloatVars.at("part_dphi")->push_back(deltaPhi(p.phi, jet->Phi));
                data.vfloatVars.at("part_d0val")->push_back(p.d0);
                data.vfloatVars.at("part_d0err")->push_back(p.d0err);
                data.vfloatVars.at("part_dzval")->push_back((pv && p.dz != 0) ? (p.dz - pv->Z) : p.dz);
                data.vfloatVars.at("part_dzerr")->push_back(p.dzerr);
                data.vintVars.at("part_charge")->push_back(p.charge);
                data.vboolVars.at("part_isElectron")->push_back(p.pid == 11 || p.pid == -11);
                data.vboolVars.at("part_isMuon")->push_back(p.pid == 13 || p.pid == -13);
                data.vboolVars.at("part_isPhoton")->push_back(p.pid == 22);
                data.vboolVars.at("part_isChargedHadron")->push_back(p.charge != 0 && !(p.pid == 11 || p.pid == -11 || p.pid == 13 || p.pid == -13));
                data.vboolVars.at("part_isNeutralHadron")->push_back(p.charge == 0 && !(p.pid == 22));
            }

            // Writing genjet features
            float min_dr = 999;
            int min_dr_index = -1;
            for (Int_t j = 0; j < branchGenJet->GetEntriesFast(); ++j) {
                const Jet *genjet = (Jet *)branchGenJet->At(j);
                float dr = deltaR(genjet, jet);

                bool is_used = std::find(genjet_used_inds.begin(), genjet_used_inds.end(), j) != genjet_used_inds.end();
                if (dr < min_dr && !is_used) {
                    min_dr = dr;
                    min_dr_index = j;
                }
            }
            if (min_dr < jetR) {
                // target genjet found
                genjet_used_inds.push_back(min_dr_index);
                const Jet *genjet = (Jet *)branchGenJet->At(min_dr_index);
                
                if (debug) {
                    std::cerr << ">> debug matched genjet (pT, eta, phi) : " << genjet->PT << "  " << genjet->Eta << "  " << genjet->Phi
                              << "   dr(genjet, jet) : " << min_dr
                              << std::endl;
                }

                data.floatVars.at("genjet_pt") = genjet->PT;
                data.floatVars.at("genjet_eta") = genjet->Eta;
                data.floatVars.at("genjet_phi") = genjet->Phi;
                data.floatVars.at("genjet_energy") = genjet->P4().Energy();
                data.floatVars.at("genjet_sdmass") = genjet->SoftDroppedP4[0].M();

                // Loop over all jet's constituents
                std::vector<ParticleInfo> genparticles;
                for (Int_t j = 0; j < genjet->Constituents.GetEntriesFast(); ++j) {
                    const TObject *object = genjet->Constituents.At(j);

                    // Check if the constituent is accessible
                    if (!object)
                        continue;

                    if (object->IsA() == GenParticle::Class()) {
                        genparticles.emplace_back((GenParticle *)object);
                    }
                    const auto &p = genparticles.back();
                    if (std::abs(p.pz) > 10000 || std::abs(p.eta) > 5 || p.pt <= 0) {
                        genparticles.pop_back();
                    }
                }
                // sort particles by pt
                std::sort(genparticles.begin(), genparticles.end(), [](const auto &a, const auto &b) { return a.pt > b.pt; });

                data.intVars["genjet_nparticles"] = genparticles.size();
                for (const auto &p : genparticles) {
                    data.vfloatVars.at("genpart_px")->push_back(p.px);
                    data.vfloatVars.at("genpart_py")->push_back(p.py);
                    data.vfloatVars.at("genpart_pz")->push_back(p.pz);
                    data.vfloatVars.at("genpart_energy")->push_back(p.energy);
                    data.vfloatVars.at("genpart_jet_deta")->push_back((jet->Eta > 0 ? 1 : -1) * (p.eta - jet->Eta));
                    data.vfloatVars.at("genpart_jet_dphi")->push_back(deltaPhi(p.phi, jet->Phi));
                    data.vfloatVars.at("genpart_z")->push_back(pv ? (p.z - pv->Z): p.z);
                    data.vfloatVars.at("genpart_t")->push_back(pv ? (p.t - pv->T): p.t);
                    float absz = std::abs(pv ? (p.z - pv->Z): p.z);
                    data.vfloatVars.at("genpart_x")->push_back(absz < 1e-10 ? 0. : p.x); // fix the precision problem for p.x and p.y
                    data.vfloatVars.at("genpart_y")->push_back(absz < 1e-10 ? 0. : p.y);
                    data.vintVars.at("genpart_pid")->push_back(p.pid);
                }
            }


            tree->Fill();
            ++num_processed;
        } // end loop of jets
    } // end loop of events

    tree->Write();
    std::cerr << TString::Format("** Written %d jets to output %s", num_processed, outputFile.Data()) << std::endl;

    delete treeReader;
    delete chain;
    delete fout;
}

//------------------------------------------------------------------------------
