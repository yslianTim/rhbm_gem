#pragma once

#include <cstdint>
#include <string>
#include <unordered_set>

#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>

class AtomSelector
{
    std::unordered_set<std::string> veto_chain_set, pick_chain_set;
    std::unordered_set<Residue> veto_residue_set, pick_residue_set;
    std::unordered_set<Element> veto_element_set, pick_element_set;

public:
    AtomSelector();
    ~AtomSelector() = default;
    std::string Describe() const;
    void Print() const;
    bool GetSelectionFlag(const std::string &, Residue, Element) const;
    void VetoChainID(const std::string & name);
    void VetoResidueType(const std::string & name);
    void VetoElementType(const std::string & name);
    void PickChainID(const std::string & name);
    void PickResidueType(const std::string & name);
    void PickElementType(const std::string & name);

};
