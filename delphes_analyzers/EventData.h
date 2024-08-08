#ifndef EventData_h
#define EventData_h

#include "TTree.h"
#include <vector>
#include <map>
#include <string>
#include <stdexcept>

// Define a struct to hold all branch variables
struct EventData {
    std::vector<std::pair<std::string, std::string>> branchList;
    std::map<std::string, bool> boolVars;
    std::map<std::string, int> intVars;
    std::map<std::string, float> floatVars;
    std::map<std::string, std::vector<bool>*> vboolVars;
    std::map<std::string, std::vector<int>*> vintVars;
    std::map<std::string, std::vector<float>*> vfloatVars;

    EventData(std::vector<std::pair<std::string, std::string>>& branchList_) {
        branchList = branchList_;
        for (const auto& pair : branchList_) {
            if (pair.second == "bool")  boolVars[pair.first] = 0;
            else if (pair.second == "int")  intVars[pair.first] = 0;
            else if (pair.second == "float")  floatVars[pair.first] = 0;
            else if (pair.second == "vector<bool>")  vboolVars[pair.first] = nullptr;
            else if (pair.second == "vector<int>")  vintVars[pair.first] = nullptr;
            else if (pair.second == "vector<float>")  vfloatVars[pair.first] = nullptr;
            else  throw std::runtime_error("Invalid branch type: " + pair.second);
        }
    }

    ~EventData() {
        for (auto& pair : vboolVars)
            if (pair.second)  delete pair.second;
        for (auto& pair : vintVars)
            if (pair.second)  delete pair.second;
        for (auto& pair : vfloatVars)
            if (pair.second)  delete pair.second;
    }

    // Reset all variables
    void reset() {
        for (auto& pair : boolVars)   pair.second = 0;
        for (auto& pair : intVars)    pair.second = 0;
        for (auto& pair : floatVars)  pair.second = 0;
        for (auto& pair : vboolVars) {
            if (pair.second)  pair.second->clear();
            else  pair.second = new std::vector<bool>;
        }
        for (auto& pair : vintVars) {
            if (pair.second)  pair.second->clear();
            else  pair.second = new std::vector<int>;
        }
        for (auto& pair : vfloatVars) {
            if (pair.second)  pair.second->clear();
            else  pair.second = new std::vector<float>;
        }
    }

    // Copy from another EventData
    void copy(const EventData& data) {
        for (auto& pair : boolVars)   pair.second = data.boolVars.at(pair.first);
        for (auto& pair : intVars)    pair.second = data.intVars.at(pair.first);
        for (auto& pair : floatVars)  pair.second = data.floatVars.at(pair.first);
        for (auto& pair : vboolVars) {
            if (pair.second)  *pair.second = *data.vboolVars.at(pair.first);
            else  vboolVars[pair.first] = new std::vector<bool>(*data.vboolVars.at(pair.first));
        }
        for (auto& pair : vintVars) {
            if (pair.second)  *pair.second = *data.vintVars.at(pair.first);
            else  vintVars[pair.first] = new std::vector<int>(*data.vintVars.at(pair.first));
        }
        for (auto& pair : vfloatVars) {
            if (pair.second)  *pair.second = *data.vfloatVars.at(pair.first);
            else  vfloatVars[pair.first] = new std::vector<float>(*data.vfloatVars.at(pair.first));
        }
    }

    // Set branch addresses for an input tree
    void setBranchAddresses(TTree* tree) {
        for (auto& pair : boolVars)     tree->SetBranchAddress(pair.first.c_str(), &pair.second);
        for (auto& pair : intVars)      tree->SetBranchAddress(pair.first.c_str(), &pair.second);
        for (auto& pair : floatVars)    tree->SetBranchAddress(pair.first.c_str(), &pair.second);
        for (auto& pair : vboolVars)    tree->SetBranchAddress(pair.first.c_str(), &pair.second);
        for (auto& pair : vintVars)     tree->SetBranchAddress(pair.first.c_str(), &pair.second);
        for (auto& pair : vfloatVars)   tree->SetBranchAddress(pair.first.c_str(), &pair.second);
    }

    // Configure branches for an output tree
    void setOutputBranch(TTree* tree) {
        // follow the branch list order
        for (const auto& pair : branchList) {
            if (pair.second == "bool") {
                tree->Branch(pair.first.c_str(), &boolVars.at(pair.first));
            }
            else if (pair.second == "int") {
                tree->Branch(pair.first.c_str(), &intVars.at(pair.first));
            }
            else if (pair.second == "float") {
                tree->Branch(pair.first.c_str(), &floatVars.at(pair.first));
            }
            else if (pair.second == "vector<bool>") {
                tree->Branch(pair.first.c_str(), &vboolVars.at(pair.first), /*bufsize=*/102400);
            }
            else if (pair.second == "vector<int>") {
                tree->Branch(pair.first.c_str(), &vintVars.at(pair.first), /*bufsize=*/102400);
            }
            else if (pair.second == "vector<float>") {
                tree->Branch(pair.first.c_str(), &vfloatVars.at(pair.first), /*bufsize=*/102400);
            }
            else {
                throw std::runtime_error("Invalid branch type: " + pair.second);
            }
        }
    }

};

#endif