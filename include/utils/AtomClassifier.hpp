#pragma once

#include <vector>
#include <tuple>
#include "AtomicInfoHelper.hpp"

class AtomClassifier
{
    using ElementKeyType = GroupKeyMapping<ElementGroupClassifierTag>::type;
    using ResidueKeyType = GroupKeyMapping<ResidueGroupClassifierTag>::type;

    static const int m_main_chain_member_count{ 4 };
    static const std::vector<Element> m_main_chain_member_element_list;
    static const std::vector<Remoteness> m_main_chain_member_remoteness_list;

public:
    AtomClassifier(void);
    ~AtomClassifier();

    static Element GetMainChainElement(size_t id);
    static Remoteness GetMainChainRemoteness(size_t id);
    ElementKeyType GetMainChainElementClassGroupKey(size_t id) const;
    std::vector<ResidueKeyType> GetMainChainResidueClassGroupKeyList(size_t id, Residue residue = Residue::UNK) const;

private:


};