#include "AtomClassifier.hpp"

using ElementKeyType = GroupKeyMapping<ElementGroupClassifierTag>::type;
using ResidueKeyType = GroupKeyMapping<ResidueGroupClassifierTag>::type;

const std::vector<Element> AtomClassifier::m_main_chain_member_element_list
{
    Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN
};

const std::vector<Remoteness> AtomClassifier::m_main_chain_member_remoteness_list
{
    Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE, Remoteness::NONE
};

AtomClassifier::AtomClassifier(void)
{

}

AtomClassifier::~AtomClassifier()
{

}

ElementKeyType AtomClassifier::GetMainChainElementClassGroupKey(size_t id)
{
    if (id >= m_main_chain_member_count)
    {
        std::cerr << "Invalid id: " << id << std::endl;
        return {};
    }
    auto element_id{ static_cast<int>(m_main_chain_member_element_list.at(id)) };
    auto remoteness_id{ static_cast<int>(m_main_chain_member_remoteness_list.at(id)) };
    return std::make_tuple(element_id, remoteness_id, false);
}

std::vector<ResidueKeyType> AtomClassifier::GetMainChainResidueClassGroupKeyList(
    size_t id, int residue)
{
    if (id >= m_main_chain_member_count)
    {
        std::cerr << "Invalid id: " << id << std::endl;
        return {};
    }
    std::vector<ResidueKeyType> group_key_list;
    auto element_id{ static_cast<int>(m_main_chain_member_element_list.at(id)) };
    auto remoteness_id{ static_cast<int>(m_main_chain_member_remoteness_list.at(id)) };

    if (residue == -1)
    {
        group_key_list.reserve(static_cast<size_t>(AtomicInfoHelper::GetStandardResidueCount()));
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