#pragma once

#include <iostream>
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

    ElementKeyType GetMainChainElementClassGroupKey(size_t id);
    std::vector<ResidueKeyType> GetMainChainResidueClassGroupKeyList(size_t id, int residue = -1);

private:


};