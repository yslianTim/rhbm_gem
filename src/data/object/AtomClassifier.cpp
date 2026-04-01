#include <rhbm_gem/data/object/AtomClassifier.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/KeyPacker.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <stdexcept>

namespace rhbm_gem {

const std::vector<Element> AtomClassifier::m_main_chain_member_element_list
{
    Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
    Element::PHOSPHORUS, Element::OXYGEN, Element::OXYGEN,
    Element::OXYGEN, Element::OXYGEN, Element::OXYGEN, Element::OXYGEN,
    Element::CARBON, Element::CARBON, Element::CARBON,
    Element::CARBON, Element::CARBON
};

const std::vector<Spot> AtomClassifier::m_main_chain_member_spot_list
{
    Spot::CA, Spot::C, Spot::N, Spot::O,
    Spot::P, Spot::OP1, Spot::OP2,
    Spot::O5p, Spot::O4p, Spot::O3p, Spot::O2p,
    Spot::C5p, Spot::C4p, Spot::C3p,
    Spot::C2p, Spot::C1p
};

bool AtomClassifier::IsMainChainMember(Spot spot, size_t & main_chain_member_id)
{
    for (size_t i = 0; i < m_main_chain_member_count; i++)
    {
        if (m_main_chain_member_spot_list.at(i) == spot)
        {
            main_chain_member_id = i;
            return true;
        }
    }
    return false;
}

size_t AtomClassifier::GetMainChainMemberCount()
{
    return m_main_chain_member_count;
}

bool AtomClassifier::IsValidMainChainMemberID(size_t member_id)
{
    if (member_id >= m_main_chain_member_count)
    {
        Logger::Log(
            LogLevel::Error,
            "Invalid main chain member ID: " + std::to_string(member_id));
        return false;
    }
    return true;
}

Element AtomClassifier::GetMainChainElement(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return Element::UNK;
    return m_main_chain_member_element_list.at(id);
}

Spot AtomClassifier::GetMainChainSpot(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return Spot::UNK;
    return m_main_chain_member_spot_list.at(id);
}

Spot AtomClassifier::GetNucleotideMainChainSpot(size_t id)
{
    id += m_nucleotide_main_chain_member_start_id;
    if (IsValidMainChainMemberID(id) == false) return Spot::UNK;
    return m_main_chain_member_spot_list.at(id);
}

GroupKey AtomClassifier::GetGroupKeyInClass(
    const AtomObject * atom_object,
    const std::string & class_key)
{
    if (class_key == ChemicalDataHelper::GetSimpleAtomClassKey())
    {
        return KeyPackerSimpleAtomClass::Pack(
            atom_object->GetAtomKey(),
            atom_object->GetSpecialAtomFlag());
    }
    if (class_key == ChemicalDataHelper::GetComponentAtomClassKey())
    {
        return KeyPackerComponentAtomClass::Pack(
            atom_object->GetComponentKey(),
            atom_object->GetAtomKey());
    }
    if (class_key == ChemicalDataHelper::GetStructureAtomClassKey())
    {
        return KeyPackerStructureAtomClass::Pack(
            atom_object->GetStructure(),
            atom_object->GetComponentKey(),
            atom_object->GetAtomKey());
    }
    throw std::runtime_error("Unsupported class key." + class_key);
}

GroupKey AtomClassifier::GetGroupKeyInClass(ComponentKey component_key, AtomKey atom_key)
{
    return KeyPackerComponentAtomClass::Pack(component_key, atom_key);
}

GroupKey AtomClassifier::GetGroupKeyInClass(
    Structure structure,
    ComponentKey component_key,
    AtomKey atom_key)
{
    return KeyPackerStructureAtomClass::Pack(structure, component_key, atom_key);
}

GroupKey AtomClassifier::GetMainChainSimpleAtomClassGroupKey(size_t id)
{
    if (IsValidMainChainMemberID(id) == false)
    {
        Logger::Log(LogLevel::Warning, "Invalid main chain member ID: " + std::to_string(id));
        return 0;
    }
    return KeyPackerSimpleAtomClass::Pack(
        static_cast<AtomKey>(m_main_chain_member_spot_list.at(id)),
        false);
}

GroupKey AtomClassifier::GetMainChainComponentAtomClassGroupKey(size_t id, Residue residue)
{
    if (IsValidMainChainMemberID(id) == false)
    {
        Logger::Log(LogLevel::Warning, "Invalid main chain member ID: " + std::to_string(id));
        return 0;
    }
    const auto component_key{ static_cast<ComponentKey>(residue) };
    return KeyPackerComponentAtomClass::Pack(
        component_key,
        static_cast<AtomKey>(m_main_chain_member_spot_list.at(id)));
}

GroupKey AtomClassifier::GetMainChainStructureAtomClassGroupKey(
    size_t id,
    Structure structure,
    Residue residue)
{
    if (IsValidMainChainMemberID(id) == false)
    {
        Logger::Log(LogLevel::Warning, "Invalid main chain member ID: " + std::to_string(id));
        return 0;
    }
    const auto component_key{ static_cast<ComponentKey>(residue) };
    return KeyPackerStructureAtomClass::Pack(
        structure,
        component_key,
        static_cast<AtomKey>(m_main_chain_member_spot_list.at(id)));
}

std::vector<GroupKey> AtomClassifier::GetMainChainComponentAtomClassGroupKeyList(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return {};
    std::vector<GroupKey> group_key_list;
    group_key_list.reserve(ChemicalDataHelper::GetStandardResidueCount());
    for (const auto residue_id : ChemicalDataHelper::GetStandardAminoAcidList())
    {
        group_key_list.emplace_back(GetMainChainComponentAtomClassGroupKey(id, residue_id));
    }
    for (const auto residue_id : ChemicalDataHelper::GetStandardNucleotideList())
    {
        group_key_list.emplace_back(GetMainChainComponentAtomClassGroupKey(id, residue_id));
    }
    return group_key_list;
}

std::vector<GroupKey> AtomClassifier::GetMainChainStructureAtomClassGroupKeyList(
    size_t id,
    Structure structure)
{
    if (IsValidMainChainMemberID(id) == false) return {};
    std::vector<GroupKey> group_key_list;
    group_key_list.reserve(ChemicalDataHelper::GetStandardResidueCount());
    for (const auto residue_id : ChemicalDataHelper::GetStandardAminoAcidList())
    {
        group_key_list.emplace_back(
            GetMainChainStructureAtomClassGroupKey(id, structure, residue_id));
    }
    for (const auto residue_id : ChemicalDataHelper::GetStandardNucleotideList())
    {
        group_key_list.emplace_back(
            GetMainChainStructureAtomClassGroupKey(id, structure, residue_id));
    }
    return group_key_list;
}

GroupKey AtomClassifier::GetNucleotideMainChainComponentAtomClassGroupKey(
    size_t id,
    Residue residue)
{
    id += m_nucleotide_main_chain_member_start_id;
    return GetMainChainComponentAtomClassGroupKey(id, residue);
}

std::vector<GroupKey> AtomClassifier::GetNucleotideMainChainComponentAtomClassGroupKeyList(
    Residue residue)
{
    std::vector<GroupKey> group_key_list;
    group_key_list.reserve(m_nucleotide_main_chain_member_count);
    for (size_t i = 0; i < m_nucleotide_main_chain_member_count; i++)
    {
        group_key_list.emplace_back(GetNucleotideMainChainComponentAtomClassGroupKey(i, residue));
    }
    return group_key_list;
}

} // namespace rhbm_gem
