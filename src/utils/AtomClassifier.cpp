#include "AtomClassifier.hpp"
#include "KeyPacker.hpp"

#include <iostream>

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

Element AtomClassifier::GetMainChainElement(size_t id)
{
    if (id >= m_main_chain_member_count)
    {
        std::cerr << "Invalid id: " << id << std::endl;
        return Element::UNK;
    }
    return m_main_chain_member_element_list.at(id);
}

Remoteness AtomClassifier::GetMainChainRemoteness(size_t id)
{
    if (id >= m_main_chain_member_count)
    {
        std::cerr << "Invalid id: " << id << std::endl;
        return Remoteness::UNK;
    }
    return m_main_chain_member_remoteness_list.at(id);
}

uint64_t AtomClassifier::GetMainChainElementClassGroupKey(size_t id) const
{
    if (id >= m_main_chain_member_count)
    {
        std::cerr << "Invalid id: " << id << std::endl;
        return {};
    }
    auto element_id{ m_main_chain_member_element_list.at(id) };
    auto remoteness_id{ m_main_chain_member_remoteness_list.at(id) };
    return KeyPackerElementClass::Pack(element_id, remoteness_id, false);
}

std::vector<uint64_t> AtomClassifier::GetMainChainResidueClassGroupKeyList(
    size_t id, Residue residue) const
{
    if (id >= m_main_chain_member_count)
    {
        std::cerr << "Invalid id: " << id << std::endl;
        return {};
    }
    std::vector<uint64_t> group_key_list;
    auto element_id{ m_main_chain_member_element_list.at(id) };
    auto remoteness_id{ m_main_chain_member_remoteness_list.at(id) };

    if (residue == Residue::UNK)
    {
        group_key_list.reserve(static_cast<size_t>(AtomicInfoHelper::GetStandardResidueCount()));
        for (auto residue_id : AtomicInfoHelper::GetStandardResidueList())
        {
            auto group_key{ KeyPackerResidueClass::Pack(residue_id, element_id, remoteness_id, Branch::NONE, false) };
            group_key_list.emplace_back(group_key);
        }
    }
    else
    {
        auto group_key{ KeyPackerResidueClass::Pack(residue, element_id, remoteness_id, Branch::NONE, false) };
        group_key_list.emplace_back(group_key);
    }

    return group_key_list;
}