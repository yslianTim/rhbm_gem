#include "ModelObject.hpp"
#include "AtomObject.hpp"
#include "DataObjectVisitorBase.hpp"
#include "GroupPotentialEntry.hpp"

ModelObject::ModelObject(void) :
    m_key_tag{ "" }, m_pdb_id{ "" }, m_emd_id{ "" }
{
}

ModelObject::ModelObject(std::vector<std::unique_ptr<AtomObject>> atom_object_list) :
    m_atom_list{ std::move(atom_object_list) },
    m_key_tag{ "" }, m_pdb_id{ "" }, m_emd_id{ "" }
{
}

ModelObject::~ModelObject()
{

}

ModelObject::ModelObject(const ModelObject & other) :
    m_key_tag{ other.m_key_tag }, m_pdb_id{ other.m_pdb_id }, m_emd_id{ other.m_emd_id }
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

int ModelObject::GetNumberOfSelectedAtom(void) const
{
    auto count{ 0 };
    for (auto & atom : m_atom_list)
    {
        if (atom->GetSelectedFlag() == true) count++;
    }
    return count;
}