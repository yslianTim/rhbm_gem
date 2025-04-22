#include "ModelObject.hpp"
#include "AtomObject.hpp"
#include "DataObjectVisitorBase.hpp"
#include "GroupPotentialEntry.hpp"
#include "ArrayStats.hpp"

ModelObject::ModelObject(void) :
    m_key_tag{ "" }, m_pdb_id{ "" }, m_emd_id{ "" }, m_kd_tree_root{ nullptr }
{
}

ModelObject::ModelObject(std::vector<std::unique_ptr<AtomObject>> atom_object_list) :
    m_atom_list{ std::move(atom_object_list) },
    m_key_tag{ "" }, m_pdb_id{ "" }, m_emd_id{ "" }, m_kd_tree_root{ nullptr }
{
}

ModelObject::~ModelObject()
{

}

ModelObject::ModelObject(const ModelObject & other) :
    m_key_tag{ other.m_key_tag }, m_pdb_id{ other.m_pdb_id }, m_emd_id{ other.m_emd_id },
    m_kd_tree_root{ nullptr }
{
    m_atom_list.clear();
    m_atom_list.reserve(other.m_atom_list.size());
    for (auto & atom : other.m_atom_list)
    {
        m_atom_list.emplace_back(atom->AtomObjectClone());
    }
}

std::unique_ptr<DataObjectBase> ModelObject::Clone() const
{
    return std::make_unique<ModelObject>(*this);
}

void ModelObject::Display(void) const
{
    std::cout << "This is ModelObject."
              << " It contains: "<< m_atom_list.size() << " atoms." << std::endl;
    std::cout << "The number of selected atom = " << GetNumberOfSelectedAtom() << std::endl;
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

void ModelObject::AddGroupPotentialEntry(const std::string & class_key, std::unique_ptr<GroupPotentialEntry> & entry)
{
    m_group_potential_entry_map[class_key] = std::move(entry);
}

size_t ModelObject::GetNumberOfSelectedAtom(void) const
{
    size_t count{ 0 };
    for (auto & atom : m_atom_list)
    {
        if (atom->GetSelectedFlag() == true) count++;
    }
    return count;
}

void ModelObject::BuildKDTreeRoot(void)
{
    if (m_kd_tree_root != nullptr) return;
    std::cout <<" ModelObject::BuildKDTreeRoot , #atom = "<< m_atom_list.size()<< std::endl;
    std::vector<AtomObject *> atom_ptr_list;
    atom_ptr_list.reserve(m_atom_list.size());
    for (auto & atom : m_atom_list)
    {
        atom_ptr_list.emplace_back(atom.get());
    }
    m_kd_tree_root = KDTreeAlgorithm<AtomObject>::BuildKDTree(atom_ptr_list, 0);
}

std::tuple<double, double> ModelObject::GetModelPositionRange(
    int axis, double margin) const
{
    std::vector<double> position_list;
    position_list.reserve(m_atom_list.size());
    for (auto & atom : m_atom_list)
    {
        position_list.emplace_back(atom->GetPosition().at(static_cast<size_t>(axis)));
    }
    return ArrayStats<double>::ComputeScalingRangeTuple(position_list, margin);
}

GroupPotentialEntry * ModelObject::GetGroupPotentialEntry(const std::string & class_key) const
{
    return m_group_potential_entry_map.at(class_key).get();
}

const std::unordered_map<std::string, std::unique_ptr<GroupPotentialEntry>> & ModelObject::GetGroupPotentialEntryMap(void) const
{
    return m_group_potential_entry_map;
}