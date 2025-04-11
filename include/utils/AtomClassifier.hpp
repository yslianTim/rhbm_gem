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
    const int m_main_chain_element_index[m_main_chain_member_count]   { 6, 6, 7, 8 };
    const int m_main_chain_remoteness_index[m_main_chain_member_count]{ 1, 0, 0, 0 };

public:
    AtomClassifier(void);
    ~AtomClassifier();

    ElementKeyType GetMainChainElementClassGroupKey(int id);
    std::vector<ResidueKeyType> GetMainChainResidueClassGroupKeyList(int id, int residue = -1);

private:


};