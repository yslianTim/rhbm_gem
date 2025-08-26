#include "ModelObject.hpp"
#include "AtomObject.hpp"
#include "DataObjectVisitorBase.hpp"
#include "GroupPotentialEntry.hpp"
#include "KDTreeAlgorithm.hpp"
#include "ArrayStats.hpp"
#include "Logger.hpp"

ModelObject::ModelObject(void) :
    m_key_tag{ "" }, m_pdb_id{ "" }, m_emd_id{ "" }, m_kd_tree_root{ nullptr }
{
}

ModelObject::ModelObject(std::vector<std::unique_ptr<AtomObject>> atom_object_list) :
    m_atom_list{ std::move(atom_object_list) },
    m_key_tag{ "" }, m_pdb_id{ "" }, m_emd_id{ "" }, m_kd_tree_root{ nullptr }
{
    Update();
}

ModelObject::~ModelObject()
{

}

ModelObject::ModelObject(const ModelObject & other) :
    m_key_tag{ other.m_key_tag }, m_pdb_id{ other.m_pdb_id }, m_emd_id{ other.m_emd_id },
    m_resolution_method{ other.m_resolution_method }, m_resolution{ other.m_resolution },
    m_kd_tree_root{ nullptr }
{
    m_atom_list.clear();
    m_atom_list.reserve(other.m_atom_list.size());
    for (auto & atom : other.m_atom_list)
    {
        m_atom_list.emplace_back(atom->AtomObjectClone());
    }
    Update();
}

std::unique_ptr<DataObjectBase> ModelObject::Clone() const
{
    return std::make_unique<ModelObject>(*this);
}

void ModelObject::Display(void) const
{
    Logger::Log(LogLevel::Info, "ModelObject Display: " + GetKeyTag());
    Logger::Log(LogLevel::Info, "This is ModelObject, it contains: "
                + std::to_string(m_atom_list.size()) + " atoms.");
}

void ModelObject::Update(void)
{
    BuildSelectedAtomList();
}

void ModelObject::Accept(DataObjectVisitorBase * visitor)
{
    for (const auto & atom : m_atom_list)
    {
        atom->Accept(visitor);
    }
    visitor->VisitModelObject(this);
}

void ModelObject::AddAtom(std::unique_ptr<AtomObject> component)
{
    m_atom_list.emplace_back(std::move(component));
}

void ModelObject::AddGroupPotentialEntry(
    const std::string & class_key, std::unique_ptr<GroupPotentialEntry> & entry)
{
    m_group_potential_entry_map[class_key] = std::move(entry);
}

void ModelObject::BuildKDTreeRoot(void)
{
    if (m_kd_tree_root != nullptr) return;
    Logger::Log(LogLevel::Debug, "Building KDTree for ModelObject..., including"
                + std::to_string(m_atom_list.size()) + " atoms.");
    std::vector<AtomObject *> atom_ptr_list;
    atom_ptr_list.reserve(m_atom_list.size());
    for (auto & atom : m_atom_list)
    {
        atom_ptr_list.emplace_back(atom.get());
    }
    m_kd_tree_root = KDTreeAlgorithm<AtomObject>::BuildKDTree(atom_ptr_list, 0);
}

std::array<float, 3> ModelObject::GetCenterOfMassPosition(void)
{
    if (m_center_of_mass_position != nullptr)
    {
        return *(m_center_of_mass_position);
    }
    std::vector<float> pos_x, pos_y, pos_z;
    pos_x.reserve(m_atom_list.size());
    pos_y.reserve(m_atom_list.size());
    pos_z.reserve(m_atom_list.size());
    for (auto & atom : m_atom_list)
    {
        const auto & pos{ atom->GetPositionRef() };
        pos_x.emplace_back(pos.at(0));
        pos_y.emplace_back(pos.at(1));
        pos_z.emplace_back(pos.at(2));
    }
    m_center_of_mass_position = std::make_unique<std::array<float, 3>>();
    m_center_of_mass_position->at(0) = ArrayStats<float>::ComputeMean(pos_x.data(), pos_x.size());
    m_center_of_mass_position->at(1) = ArrayStats<float>::ComputeMean(pos_y.data(), pos_y.size());
    m_center_of_mass_position->at(2) = ArrayStats<float>::ComputeMean(pos_z.data(), pos_z.size());
    return *(m_center_of_mass_position);
}

std::tuple<double, double> ModelObject::GetModelPositionRange(int axis)
{
    if (m_model_position_range[axis] != nullptr)
    {
        return *(m_model_position_range[axis]);
    }
    std::vector<double> position_list;
    position_list.reserve(m_atom_list.size());
    for (auto & atom : m_atom_list)
    {
        position_list.emplace_back(atom->GetPosition().at(static_cast<size_t>(axis)));
    }
    m_model_position_range[axis] = std::make_unique<std::tuple<double, double>>(
        ArrayStats<double>::ComputeScalingRangeTuple(position_list, 0.0));
    return *(m_model_position_range[axis]);
}

double ModelObject::GetModelPosition(int axis, double normalized_pos)
{
    auto range_tuple{ GetModelPositionRange(axis) };
    auto pos_min{ std::get<0>(range_tuple) };
    auto pos_max{ std::get<1>(range_tuple) };
    return pos_min + normalized_pos * (pos_max - pos_min);
}

double ModelObject::GetModelLength(int axis)
{
    auto range_tuple{ GetModelPositionRange(axis) };
    return std::get<1>(range_tuple) - std::get<0>(range_tuple);
}

GroupPotentialEntry * ModelObject::GetGroupPotentialEntry(const std::string & class_key) const
{
    return m_group_potential_entry_map.at(class_key).get();
}

const std::unordered_map<std::string, std::unique_ptr<GroupPotentialEntry>> &
ModelObject::GetGroupPotentialEntryMap(void) const
{
    return m_group_potential_entry_map;
}

void ModelObject::BuildSelectedAtomList(void)
{
    m_selected_atom_list.clear();
    m_selected_atom_list.reserve(m_atom_list.size());
    for (auto & atom : m_atom_list)
    {
        if (atom->GetSelectedFlag() == false) continue;
        m_selected_atom_list.emplace_back(atom.get());
    }
}

void ModelObject::FilterAtomFromSymmetry(bool is_asymmetry)
{
    if (is_asymmetry == true)
    {
        return;
    }

    for (auto & atom : m_atom_list)
    {
        auto original_selection_flag{ atom->GetSelectedFlag() };
        auto chain_id{ atom->GetChainID() };
        bool in_candidate_chain{ false };
        for (auto & [entity_id, chain_id_list] : m_chain_id_list_map)
        {
            if (chain_id == chain_id_list.at(0))
            {
                atom->SetSelectedFlag(original_selection_flag);
                in_candidate_chain = true;
                break;
            }
        }
        if (in_candidate_chain == false) atom->SetSelectedFlag(false);
    }
}
