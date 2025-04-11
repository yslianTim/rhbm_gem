#include "AtomClassifier.hpp"

using ElementKeyType = GroupKeyMapping<ElementGroupClassifierTag>::type;
using ResidueKeyType = GroupKeyMapping<ResidueGroupClassifierTag>::type;

AtomClassifier::AtomClassifier(void)
{

}

AtomClassifier::~AtomClassifier()
{

}

ElementKeyType AtomClassifier::GetMainChainElementClassGroupKey(int id)
{
    if (id >= m_main_chain_member_count)
    {
        std::cerr << "Invalid id: " << id << std::endl;
        return {};
    }
    return std::make_tuple(m_main_chain_element_index[id], m_main_chain_remoteness_index[id], false);
}

std::vector<ResidueKeyType> AtomClassifier::GetMainChainResidueClassGroupKeyList(
    int id, int residue)
{
    if (id >= m_main_chain_member_count)
    {
        std::cerr << "Invalid id: " << id << std::endl;
        return {};
    }
    std::vector<ResidueKeyType> group_key_list;
    auto element_id{ m_main_chain_element_index[id] };
    auto remoteness_id{ m_main_chain_remoteness_index[id] };

    if (residue == -1)
    {
        group_key_list.reserve(AtomicInfoHelper::GetStandardResidueCount());
        for (auto residue_id : AtomicInfoHelper::GetStandardResidueList())
        {
            auto group_key{ std::make_tuple(residue_id, element_id, remoteness_id, 0, false) };
            group_key_list.emplace_back(group_key);
        }
    }
    else
    {
        auto group_key{ std::make_tuple(residue, element_id, remoteness_id, 0, false) };
        group_key_list.emplace_back(group_key);
    }

    return group_key_list;
}