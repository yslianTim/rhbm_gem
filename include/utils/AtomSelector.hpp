#pragma once

#include <iostream>
#include <string>
#include <set>

class AtomSelector
{
    std::set<std::string> veto_chain_set, pick_chain_set;
    std::set<std::string> veto_indicator_set, pick_indicator_set;
    std::set<int> veto_residue_set, pick_residue_set;
    std::set<int> veto_element_set, pick_element_set;
    std::set<int> veto_remoteness_set, pick_remoteness_set;
    std::set<int> veto_branch_set, pick_branch_set;

public:
    AtomSelector(void);
    ~AtomSelector();
    void Print(void) const;
    bool GetSelectionFlag(const std::string &, const std::string &, int, int, int, int);
    void VetoChainID(const std::string & name);
    void VetoIndicator(const std::string & name);
    void VetoResidueType(const std::string & name);
    void VetoElementType(const std::string & name);
    void VetoRemotenessType(const std::string & name);
    void VetoBranchType(const std::string & name);
    void PickChainID(const std::string & name);
    void PickIndicator(const std::string & name);
    void PickResidueType(const std::string & name);
    void PickElementType(const std::string & name);
    void PickRemotenessType(const std::string & name);
    void PickBranchType(const std::string & name);

};