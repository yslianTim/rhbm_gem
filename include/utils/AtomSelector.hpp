#pragma once

#include <cstdint>
#include <string>
#include <unordered_set>

enum class Residue : uint16_t;
enum class Element : uint16_t;
enum class Remoteness : uint8_t;

class AtomSelector
{
    std::unordered_set<std::string> veto_chain_set, pick_chain_set;
    std::unordered_set<Residue> veto_residue_set, pick_residue_set;
    std::unordered_set<Element> veto_element_set, pick_element_set;
    std::unordered_set<Remoteness> veto_remoteness_set, pick_remoteness_set;

public:
    AtomSelector(void);
    ~AtomSelector();
    void Print(void) const;
    bool GetSelectionFlag(const std::string &, Residue, Element, Remoteness) const;
    void VetoChainID(const std::string & name);
    void VetoResidueType(const std::string & name);
    void VetoElementType(const std::string & name);
    void VetoRemotenessType(const std::string & name);
    void PickChainID(const std::string & name);
    void PickResidueType(const std::string & name);
    void PickElementType(const std::string & name);
    void PickRemotenessType(const std::string & name);

};
