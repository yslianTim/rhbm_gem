#include "internal/io/file/CifFormat.hpp"

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/AtomicModelDataBlock.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/ChemicalComponentEntry.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/ComponentHelper.hpp>
#include <rhbm_gem/utils/math/KDTreeAlgorithm.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <cstdint>
#include <unordered_set>

namespace rhbm_gem {

namespace
{

struct CanonicalAtomPair
{
    const AtomObject * atom_a;
    const AtomObject * atom_b;

    bool operator==(const CanonicalAtomPair & other) const
    {
        return atom_a == other.atom_a && atom_b == other.atom_b;
    }
};

struct CanonicalAtomPairHash
{
    size_t operator()(const CanonicalAtomPair & pair) const
    {
        auto a{ reinterpret_cast<std::uintptr_t>(pair.atom_a) };
        auto b{ reinterpret_cast<std::uintptr_t>(pair.atom_b) };
        return std::hash<std::uintptr_t>{}(a ^ (b << 1));
    }
};

CanonicalAtomPair BuildCanonicalAtomPair(const AtomObject * atom_1, const AtomObject * atom_2)
{
    if (atom_1 < atom_2) return {atom_1, atom_2};
    return {atom_2, atom_1};
}

} // namespace

void CifFormat::ConstructBondList()
{
    BuildPepetideBondEntry();
    BuildPhosphodiesterBondEntry();
    const auto bond_count_before{ m_data_block->GetBondObjectList().size() };
    const auto & atom_object_map{ m_data_block->GetAtomObjectMap() };
    auto bond_key_system{ m_data_block->GetBondKeySystemPtr() };
    auto component_key_system{ m_data_block->GetComponentKeySystemPtr() };
    for (const auto & [model_number, atom_object_list] : atom_object_map)
    {
        if (atom_object_list.empty()) continue;
        std::vector<AtomObject *> atom_ptr_list;
        atom_ptr_list.reserve(atom_object_list.size());
        for (const auto & atom : atom_object_list)
        {
            atom_ptr_list.emplace_back(atom.get());
        }
        auto kd_tree_root{ KDTreeAlgorithm<AtomObject>::BuildKDTree(atom_ptr_list, 0) };
        std::unordered_set<CanonicalAtomPair, CanonicalAtomPairHash> processed_bond_pair_set;
        processed_bond_pair_set.reserve(atom_object_list.size() * 2);
        for (const auto & atom : atom_object_list)
        {
            auto component_id_1{ atom->GetComponentID() };
            auto atom_id_1{ atom->GetAtomID() };
            auto sequence_id_1{ atom->GetSequenceID() };
            auto chain_id_1{ atom->GetChainID() };
            auto neighbor_atom_list{
                KDTreeAlgorithm<AtomObject>::RangeSearch(
                    kd_tree_root.get(), atom.get(), m_bond_searching_radius)
            };

            for (auto neighbor_atom : neighbor_atom_list)
            {
                if (neighbor_atom == atom.get()) continue;
                auto canonical_pair{ BuildCanonicalAtomPair(atom.get(), neighbor_atom) };
                if (processed_bond_pair_set.find(canonical_pair) != processed_bond_pair_set.end())
                {
                    continue;
                }
                processed_bond_pair_set.insert(canonical_pair);

                auto component_id_2{ neighbor_atom->GetComponentID() };
                auto atom_id_2{ neighbor_atom->GetAtomID() };
                auto sequence_id_2{ neighbor_atom->GetSequenceID() };
                auto chain_id_2{ neighbor_atom->GetChainID() };
                if (bond_key_system->IsRegistedBond(atom_id_1, atom_id_2) == false) continue;
                auto component_key_1{ component_key_system->GetComponentKey(component_id_1) };
                auto bond_key{ bond_key_system->GetBondKey(atom_id_1, atom_id_2) };
                if (m_data_block->HasComponentBondEntry(component_key_1, bond_key) == false) continue;

                bool is_in_same_component{ (component_id_1 == component_id_2) };
                bool is_in_same_chain{ (chain_id_1 == chain_id_2) };
                bool is_in_consecutive_sequence{ (sequence_id_1 + 1 == sequence_id_2) };

                // Peptide bond C-N between consecutive residues in the same chain
                bool is_peptide_bond{
                    !is_in_same_component &&
                    is_in_same_chain &&
                    is_in_consecutive_sequence &&
                    (bond_key == static_cast<BondKey>(Link::C_N))
                };

                // Phosphodiester bond P-O3' between consecutive residues in the same chain
                bool is_phosphodiester_bond{
                    !is_in_same_component &&
                    is_in_same_chain &&
                    is_in_consecutive_sequence &&
                    (bond_key == static_cast<BondKey>(Link::P_O3p))
                };

                if (is_in_same_component == false &&
                    is_peptide_bond == false &&
                    is_phosphodiester_bond == false) continue;

                auto bond_entry{
                    m_data_block->GetComponentBondEntryPtr(component_key_1, bond_key)
                };
                if (bond_entry == nullptr) continue;

                auto bond_object{ std::make_unique<BondObject>(atom.get(), neighbor_atom) };
                bond_object->SetBondKey(bond_key);
                bond_object->SetBondType(bond_entry->bond_type);
                bond_object->SetBondOrder(bond_entry->bond_order);
                bond_object->SetSpecialBondFlag(false);
                m_data_block->AddBondObject(std::move(bond_object));
            }
        }
    }
    Logger::Log(LogLevel::Info,
        "Construct " + std::to_string(m_data_block->GetBondObjectList().size() - bond_count_before)
        + " bonds (total " + std::to_string(m_data_block->GetBondObjectList().size()) + ").");
}

void CifFormat::BuildDefaultChemicalComponentEntry(const std::string & comp_id)
{
    auto entry{ std::make_unique<ChemicalComponentEntry>() };
    entry->SetComponentId(comp_id);
    entry->SetComponentName("");
    entry->SetComponentType("");
    entry->SetComponentFormula("");
    entry->SetComponentMolecularWeight(0.0f);
    auto standard_flag{ false };
    if (ChemicalDataHelper::GetResidueFromString(comp_id) != Residue::UNK) standard_flag = true;
    entry->SetStandardMonomerFlag(standard_flag);

    m_data_block->GetComponentKeySystemPtr()->RegisterComponent(comp_id);
    auto component_key{ m_data_block->GetComponentKeySystemPtr()->GetComponentKey(comp_id) };
    m_data_block->AddChemicalComponentEntry(component_key, std::move(entry));
}

void CifFormat::BuildDefaultComponentAtomEntry(
    const std::string & comp_id,
    const std::string & atom_id,
    const std::string & element_symbol)
{
    auto component_key{ m_data_block->GetComponentKeySystemPtr()->GetComponentKey(comp_id) };

    m_data_block->GetAtomKeySystemPtr()->RegisterAtom(atom_id);
    auto atom_key{ m_data_block->GetAtomKeySystemPtr()->GetAtomKey(atom_id) };

    ComponentAtomEntry atom_entry;
    atom_entry.atom_id = atom_id;
    atom_entry.element_type = ChemicalDataHelper::GetElementFromString(element_symbol);
    atom_entry.aromatic_atom_flag = false;
    atom_entry.stereo_config = StereoChemistry::NONE;

    m_data_block->AddComponentAtomEntry(component_key, atom_key, atom_entry);
}

void CifFormat::BuildDefaultComponentBondEntry()
{
    for (auto & residue : ChemicalDataHelper::GetStandardAminoAcidList())
    {
        auto comp_id{ ChemicalDataHelper::GetLabel(residue) };
        auto component_key{ m_data_block->GetComponentKeySystemPtr()->GetComponentKey(comp_id) };
        for (auto & link : ComponentHelper::GetLinkList(residue))
        {
            auto bond_id{ ChemicalDataHelper::GetLabel(link) };
            m_data_block->GetBondKeySystemPtr()->RegisterBond(bond_id);
            auto bond_key{ m_data_block->GetBondKeySystemPtr()->GetBondKey(bond_id) };

            ComponentBondEntry bond_entry;
            bond_entry.bond_id = bond_id;
            bond_entry.bond_type = BondType::COVALENT;
            bond_entry.bond_order = BondOrder::UNK;
            bond_entry.aromatic_atom_flag = false;
            bond_entry.stereo_config = StereoChemistry::NONE;

            m_data_block->AddComponentBondEntry(component_key, bond_key, bond_entry);
        }
    }
}

void CifFormat::BuildPepetideBondEntry()
{
    for (auto & residue : ChemicalDataHelper::GetStandardAminoAcidList())
    {
        auto comp_id{ ChemicalDataHelper::GetLabel(residue) };
        auto component_key{ m_data_block->GetComponentKeySystemPtr()->GetComponentKey(comp_id) };
        if (m_data_block->HasChemicalComponentEntry(component_key) == false) continue;
        auto bond_id{ ChemicalDataHelper::GetLabel(Link::C_N) };
        m_data_block->GetBondKeySystemPtr()->RegisterBond(bond_id);
        auto bond_key{ m_data_block->GetBondKeySystemPtr()->GetBondKey(bond_id) };

        ComponentBondEntry bond_entry;
        bond_entry.bond_id = bond_id;
        bond_entry.bond_type = BondType::COVALENT;
        bond_entry.bond_order = BondOrder::SINGLE;
        bond_entry.aromatic_atom_flag = false;
        bond_entry.stereo_config = StereoChemistry::NONE;

        m_data_block->AddComponentBondEntry(component_key, bond_key, bond_entry);
    }
}

void CifFormat::BuildPhosphodiesterBondEntry()
{
    for (auto & residue : ChemicalDataHelper::GetStandardNucleotideList())
    {
        auto comp_id{ ChemicalDataHelper::GetLabel(residue) };
        auto component_key{ m_data_block->GetComponentKeySystemPtr()->GetComponentKey(comp_id) };
        if (m_data_block->HasChemicalComponentEntry(component_key) == false) continue;
        auto bond_id{ ChemicalDataHelper::GetLabel(Link::P_O3p) };
        m_data_block->GetBondKeySystemPtr()->RegisterBond(bond_id);
        auto bond_key{ m_data_block->GetBondKeySystemPtr()->GetBondKey(bond_id) };

        ComponentBondEntry bond_entry;
        bond_entry.bond_id = bond_id;
        bond_entry.bond_type = BondType::COVALENT;
        bond_entry.bond_order = BondOrder::SINGLE;
        bond_entry.aromatic_atom_flag = false;
        bond_entry.stereo_config = StereoChemistry::NONE;

        m_data_block->AddComponentBondEntry(component_key, bond_key, bond_entry);
    }
}

} // namespace rhbm_gem
